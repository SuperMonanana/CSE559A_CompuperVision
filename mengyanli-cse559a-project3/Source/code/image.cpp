
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: image.cpp                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"



Image::Image()
:
Vector(),
width(0),
height(0),
colors(3)
{
	// empty
}

Image::Image(int colors)
:
Vector(),
width(0),
height(0),
colors(colors)
{
	// empty
}

Image::Image(int height, int width)
:
Vector(width*height*3),
width(width),
height(height),
colors(3)
{
	// empty
}

Image::Image(int width, int height, int colors)
:
Vector(width*height*colors),
width(width),
height(height),
colors(colors)
{
	//empty
}

int Image::getWidth() const
{
	return width;
}

int Image::getHeight() const
{
	return height;
}

int Image::getColors() const
{
	return colors;
}

void Image::resize(int new_width, int new_height)
{
	Vector::resize(new_width * new_height * colors);
	width=new_width;
	height=new_height;
}

void Image::resize(int new_width, int new_height, int new_colors)
{
	Vector::resize(new_width * new_height * new_colors);
	width=new_width;
	height=new_height;
	colors=new_colors;
}

void Image::crop(int x_min, int y_min, int x_max, int y_max, Image& result) const
{
	result.resize(x_max-x_min+1, y_max-y_min+1, result.colors);
	for (int x=x_min; x<=x_max; x++) {
		for (int y=y_min; y<=y_max; y++) {
			for (int c=0; c<colors; c++) {
				result.pixel(x-x_min, y-y_min, c) = pixel(x, y, c);
			}
		}
	}
}

void Image::line(int x_min, int y_min, int x_max, int y_max, double r, double g, double b)
{
	if (x_min==x_max) {
		for (int y=min(y_min, y_max); y<=max(y_max, y_min); y++) {
			pixel(x_min, y, 0) = b;
			pixel(x_min, y, 1) = g;
			pixel(x_min, y, 2) = b;
		}
	}
	else if (y_min==y_max) {
		for (int x=min(x_min, x_max); x<=max(x_max, x_min); x++) {
			pixel(x, y_min, 0) = b;
			pixel(x, y_min, 1) = g;
			pixel(x, y_min, 2) = b;
		}
	}
	else {
		throw Error("Line: This routine only works for horizontal and vertical lines so far");
	}
}

void Image::grayscale(Image& dest) const
{
  if (&dest == this)
    throw Error("grayscale: Source and destination image must be different.");
	dest.resize(width, height, 1);
	for (int y=0; y<height; y++) {
		for (int x=0; x<width; x++) {
			double sum=0.0;
			for (int c=0; c<colors; c++) {
				sum += pixel(x,y,c)*pixel(x,y,c);
			}
			dest.pixel(x, y, 0) = sqrt(sum/3.0);
		}
	}
}

void Image::resample(Image& dest) const
{
  if (&dest == this)
    throw Error("resample: Source and destination image must be different.");
	if (dest.colors != colors) {
		throw Error("Resample: Image destination and source must have same color depths");
	}
	double width_ratio = double(width) / double(dest.width);
	double height_ratio = double(height) / double(dest.height);
	if (width_ratio>1.3 || height_ratio>1.3) {
		double new_width_ratio;
		double new_height_ratio;
		if (width_ratio>height_ratio) {
			new_width_ratio=1.3;
			new_height_ratio=height_ratio * (new_width_ratio/width_ratio);
		}
		else {
			new_height_ratio=1.3;
			new_width_ratio=width_ratio * (new_height_ratio/height_ratio);
		}
		Image new_dest(int(double(width)/new_width_ratio+1.0), int(double(height)/new_height_ratio+1.0), colors);
		resample(new_dest);
		new_dest.resample(dest);
		return;
	}
	for (int y=0; y<dest.height; y++) {
		double source_y = height_ratio * double(y);
		for (int x=0; x<dest.width; x++) {
			double source_x = width_ratio * double(x);
			for (int c=0; c<colors; c++) {
				dest.pixel(x, y, c) = sample(source_x, source_y, c);
			}
		}
	}
}




double Image::pixelFilter(const FilterKernel& kernel, int x, int y, int c) const
{
	double result=0.0;
	for (int b=0; b<kernel.getHeight(); b++) {
		int j=y+b-kernel.getHeight()/2;
		if (j>=height) j=height-j%height-1;
		for (int a=0; a<kernel.getWidth(); a++) {
			int i=x+a-kernel.getWidth()/2;
			if (i>=width) i=width-i%width-1;
			result += kernel.pixel(a, b) * pixel(abs(i), abs(j), c);
		}
	}
	return result;
}

void Image::imageFilter(Image& dest, const FilterKernel& kernel) const
{
  if (&dest == this)
    throw Error("imageFilter: Source and destination image must be different.");
	dest.resize(width, height, colors);
	for (int y=0; y<height; y++) {
		for (int x=0; x<width; x++) {
			for (int c=0; c<colors; c++) {
				dest.pixel(x, y, c) = pixelFilter(kernel, x, y, c);
			}
		}
	}
}


void Image::add(const Image& rhs, Image& result) const
{
	if (rhs.width!=width || rhs.height!=height || rhs.colors!=colors) {
		throw Error("add: image dimentions must be the same");
	}
	result.resize(width, height, colors);
	Vector::add(rhs, result);
}

void Image::sub(const Image& rhs, Image& result) const
{
	if (rhs.width!=width || rhs.height!=height || rhs.colors!=colors) {
		throw Error("sub: image dimentions must be the same");
	}
	result.resize(width, height, colors);
	Vector::sub(rhs, result);
}



void Image::normalize(double low, double high, Image& result) const
{

	result.resize(width, height, colors);
	double range = high-low;
	double largest=-DBL_MAX;
	double smallest=DBL_MAX;
	for (int y=0; y<height; y++) {
		for (int x=0; x<width; x++) {
			for (int c=0; c<colors; c++) {
				if (pixel(x, y, c)>largest) largest = pixel(x, y, c);
				if (pixel(x, y, c)<smallest) smallest = pixel(x, y, c);
			}
		}
	}
	double distance = largest-smallest;
	if (distance <= DBL_MIN && distance >= -DBL_MIN) distance=1.0; 
	for (int y=0; y<height; y++) {
		for (int x=0; x<width; x++) {
			for (int c=0; c<colors; c++) {
				result.pixel(x, y, c) = (pixel(x, y, c) - smallest)/distance*range + low;
			}
		}
	}
}

void Image::saveTarga(std::string filename) const
{
	BinaryFileWriter file(filename);
	Targa targa;
	save(file, targa);
}

void Image::saveTarga(BinaryFileWriter& file) const
{
	Targa targa;
	save(file, targa);
}

void Image::save(std::string filename, ImageFormat& format) const
{
	BinaryFileWriter file(filename);
	save(file, format);
}

void Image::save(BinaryFileWriter& file, ImageFormat& format) const
{
	format.save(file, *this);
}

void Image::loadTarga(std::string filename)
{
	BinaryFileReader file(filename);
	Targa targa;
	load(file, targa);
}

void Image::loadTarga(BinaryFileReader& file)
{
	Targa targa;
	load(file, targa);
}

void Image::load(std::string filename, ImageFormat& format)
{
	BinaryFileReader file(filename);
	load(file, format);
}

void Image::load(BinaryFileReader& file, ImageFormat& format)
{
	format.load(file, *this);
}






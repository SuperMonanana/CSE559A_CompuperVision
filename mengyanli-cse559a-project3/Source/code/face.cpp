
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: face.cpp                                                                             //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"




Face::Face(int width, int height)
:
Image(width, height, 1)
{
}

Face::Face()
:
Image(1)
{
}

void Face::subimage(int x_low, int x_high, int y_low, int y_high, const Image& img, bool debug)
{
	static Image temp;
	temp.resize(getWidth(), getHeight(), img.getColors());
	if (x_high-x_low+1 != getWidth() || y_high-y_low+1 != getHeight()) {
		throw Error("Programmer error: subimage window is not the same size as face");
	}
	for (int y=y_low; y<=y_high; y++) {
		for (int x=x_low; x<=x_high; x++) {
			for (int c=0; c<img.getColors(); c++) {
				temp.pixel(x-x_low, y-y_low, c) = img.pixel(x, y, c);
			}
		}
	}
	static Image grayscale;
	temp.grayscale(grayscale);
	grayscale.normalize(0.0, 255.0, *this);
}

void Face::load(BinaryFileReader &file)
{
	int new_width=file.readInt();
	int new_height=file.readInt();
	resize(new_width, new_height, 1);
	int size=width*height;
	for (int i=0; i<size; i++) {
		(*this)[i] = file.readDouble();
	}
}

void Face::save(BinaryFileWriter &file) const
{
	file.write(width);
	file.write(height);
	int size=width*height;
	for (int i=0; i<size; i++) {
		file.write((*this)[i]);
	}
}



void Face::load(BinaryFileReader& file, ImageFormat& format)
{
	if (width<=0 || height<=0) {
		throw Error("Face width and height must be set before loading");
	}
	Image full_image;
	full_image.load(file, format);
	Image grayscale;
	full_image.grayscale(grayscale);
	grayscale.resample(*this);
	(*this).normalize(0.0, 255.0, (*this));
}

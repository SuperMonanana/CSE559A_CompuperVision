
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: targa.cpp                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"



Targa::Targa()
{
	//empty
}

void Targa::load(FileReader& file, Image& img)
{
	// Number of Characters in Identification Field.
        file.read(identificationFieldLength);
	// Color Map Type 0: none included 1: included
	file.read(colorMapType);
	// 2 for DATA TYPE 2: Unmapped RGB
	file.read(imageTypeCode);
//	if (imageTypeCode!=2) {
//		throw Error("Load Targa: Image Type Code is %%, expected 2", imageTypeCode);
//	}


	// *** Color Map Specification ***
	// Color Map Origin, index of first color map entry.
	file.read(colorMapOrigin);
	// Color Map Length, count of color map entries
	file.read(colorMapLength);
	// Color Map Entry Size, number of bits in color map entry
	file.read(bitsPerColorMapEntry);
	
	// *** Image Specification ***
	// X Origin of Image, X coordinate of the lower left corner of the image
	file.read(xOrigin);
	// Y Origin of Image, Y coordinate of the lower left corner of the image
	file.read(yOrigin);
	// Image Width in pixels
	file.read(width);
	//std::cout << "Image width: " << width << std::endl;
	// Image Height in pixels
	file.read(height);
	//std::cout << "Image height: " << height << std::endl;
	// Image Pixel size (bits per pixel)
	file.read(bitsPerPixel);
	if (bitsPerPixel!=32 && bitsPerPixel!=24) {
		throw Error("Load Targa: Image is %% bits per pixel, expected 24 or 32", bitsPerPixel);
	}
	// Image Descriptor Byte
	file.read(imageDescriptorByte);
//	if (imageDescriptorByte!=bitsPerPixel-24) {
//		throw Error("Load Targa: Image descriptor is %%, expected %%", imageDescriptorByte, bitsPerPixel-24);
//	}

	// skip over identification field
	file.skip(identificationFieldLength);

	// skip over color map data
	file.skip(colorMapLength*bitsPerColorMapEntry/8);

	img.resize(width, height, 3);
	// read the image data
	for (int y=height-1; y>=0; y--) {
		for (int x=0; x<width; x++) {
			for (int c=0; c<3; c++) {
				unsigned char pixel;
				file.read(pixel);
				img.pixel(x, y, 2-c) = double(pixel);
			}
			if (bitsPerPixel==32) {
				unsigned char attribute;
				file.read(attribute); // read and ignore the attribute
			}
		}
	}



}

void Targa::save(FileWriter& file, const Image& img)
{
	// Number of Characters in Identification Field.
	identificationFieldLength=0;
	file.write(identificationFieldLength);
	// Color Map Type 0: none included 1: included
	colorMapType=0;
	file.write(colorMapType);
	// 2 for DATA TYPE 2: Unmapped RGB
	imageTypeCode=2;
	file.write(imageTypeCode);

	// *** Color Map Specification ***
	// Color Map Origin, index of first color map entry.
	colorMapOrigin=0;
	file.write(colorMapOrigin);
	// Color Map Length, count of color map entries
	colorMapLength=0;
	file.write(colorMapLength);
	// Color Map Entry Size, number of bits in color map entry
	bitsPerColorMapEntry=0;
	file.write(bitsPerColorMapEntry);
	
	// *** Image Specification ***
	// X Origin of Image, X coordinate of the lower left corner of the image
	xOrigin = 0;
	file.write(xOrigin);
	// Y Origin of Image, Y coordinate of the lower left corner of the image
	yOrigin = 0;
	file.write(yOrigin);
	// Image Width in pixels
	width=img.getWidth();
	file.write(width);
	// Image Height in pixels
	height=img.getHeight();
	file.write(height);
	// Image Pixel size (bits per pixel)
	bitsPerPixel=24;
	file.write(bitsPerPixel);
	// Image Descriptor Byte
	imageDescriptorByte = 0;
	file.write(imageDescriptorByte);

	// skip identification field

	// skip color map data

	// write the image data
	for (int y=height-1; y>=0; y--) {
		for (int x=0; x<width; x++) {
			for (int c=0; c<3; c++) {
				unsigned char pixel;
				double double_value;
				if (img.getColors()==1) {
					double_value = img.pixel(x, y, 0);
				}
				else {
					double_value = img.pixel(x, y, 2-c);
				}
				int img_value = int(double_value);
				pixel = (unsigned char) (min (max(0, img_value), 255));
				file.write(pixel);
			}

		}
	}
	std::cout << "Saved Targa image '" << file.getFilename() << "'" << std::endl;


}



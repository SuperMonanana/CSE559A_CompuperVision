
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: targa.h                                                                              //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class Targa : public ImageFormat
{
public:
	Targa();
	virtual void load(FileReader& file, Image& img);
	virtual void save(FileWriter& file, const Image& img);

private:
	// Number of Characters in Identification Field.
	unsigned char identificationFieldLength;
	// Color Map Type 0: none included 1: included
	unsigned char colorMapType;
	// 2 for DATA TYPE 2: Unmapped RGB
	unsigned char imageTypeCode;


	// *** Color Map Specification ***
	// Color Map Origin, index of first color map entry.
	unsigned short colorMapOrigin;
	// Color Map Length, count of color map entries
	unsigned short colorMapLength;
	// Color Map Entry Size, number of bits in color map entry
	unsigned char bitsPerColorMapEntry;
	
	// *** Image Specification ***
	// X Origin of Image, X coordinate of the lower left corner of the image
	unsigned short xOrigin;
	// Y Origin of Image, Y coordinate of the lower left corner of the image
	unsigned short yOrigin;
	// Image Width in pixels
	unsigned short width;
	// Image Height in pixels
	unsigned short height;
	// Image Pixel size (bits per pixel)
	unsigned char bitsPerPixel;
	// Image Descriptor Byte
	unsigned char imageDescriptorByte;
	

};

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: imageformat.cpp                                                                      //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"



void ImageFormat::load(std::string filename, Image& img)
{
	BinaryFileReader file(filename);
	load(file, img);
}

void ImageFormat::save(std::string filename, const Image& img)
{
	BinaryFileWriter file(filename);
	save(file, img);
}
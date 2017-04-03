
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: filterkernel.cpp                                                                     //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"



// FilterKernel stuff
FilterKernel::FilterKernel(int height, int width)
:
Image(height, width, 1)
{
	//empty
}

double& FilterKernel::pixel(int i, int j)
{
	return Image::pixel(i, j, 0);
}

double FilterKernel::pixel(int x, int y) const
{
	return Image::pixel(x, y, 0);
}

int FilterKernel::getHeight() const
{
	return Image::getHeight();
}

int FilterKernel::getWidth() const
{
	return Image::getWidth();
}

// FilterKernel3x3 stuff
FilterKernel3x3::FilterKernel3x3()
:
FilterKernel(3, 3)
{
	//empty
}

// FilterKernel5x5 stuff
FilterKernel5x5::FilterKernel5x5()
:
FilterKernel(5, 5)
{
	//empty
}

// BlurOperator3x3 stuff

BlurOperator3x3 BlurOperator3x3::theInstance;

BlurOperator3x3::BlurOperator3x3()
:
FilterKernel3x3()
{
	// left hand column
	pixel(0, 0) = 1.0/16.0;
	pixel(0, 1) = 2.0/16.0;
	pixel(0, 2) = 1.0/16.0;

	// middle column
	pixel(1, 0) = 2.0/16.0;
	pixel(1, 1) = 4.0/16.0;
	pixel(1, 2) = 2.0/16.0;

	// right hand column
	pixel(2, 0) = 1.0/16.0;
	pixel(2, 1) = 2.0/16.0;
	pixel(2, 2) = 1.0/16.0;
}

// BlurOperator5x5 stuff

BlurOperator5x5 BlurOperator5x5::theInstance;

BlurOperator5x5::BlurOperator5x5()
:
FilterKernel5x5()
{
	// far left hand column
	pixel(0, 0) = 1.0/44.0;
	pixel(0, 1) = 1.0/44.0;
	pixel(0, 2) = 2.0/44.0;
	pixel(0, 3) = 1.0/44.0;
	pixel(0, 4) = 1.0/44.0;

	// near left hand column
	pixel(1, 0) = 1.0/44.0;
	pixel(1, 1) = 2.0/44.0;
	pixel(1, 2) = 3.0/44.0;
	pixel(1, 3) = 2.0/44.0;
	pixel(1, 4) = 1.0/44.0;

	// middle column
	pixel(2, 0) = 2.0/44.0;
	pixel(2, 1) = 3.0/44.0;
	pixel(2, 2) = 4.0/44.0;
	pixel(2, 3) = 3.0/44.0;
	pixel(2, 4) = 2.0/44.0;

	// near right hand column
	pixel(3, 0) = 1.0/44.0;
	pixel(3, 1) = 2.0/44.0;
	pixel(3, 2) = 3.0/44.0;
	pixel(3, 3) = 2.0/44.0;
	pixel(3, 4) = 1.0/44.0;

	// far right hand column
	pixel(4, 0) = 1.0/44.0;
	pixel(4, 1) = 1.0/44.0;
	pixel(4, 2) = 2.0/44.0;
	pixel(4, 3) = 1.0/44.0;
	pixel(4, 4) = 1.0/44.0;
}

const BlurOperator5x5& BlurOperator5x5::instance()
{
	return theInstance;
}

// SobelOperatorX stuff

SobelOperatorX SobelOperatorX::theInstance;

SobelOperatorX::SobelOperatorX()
:
FilterKernel3x3()
{
	// left hand column
	pixel(0, 0) = -1.0/8.0;
	pixel(0, 1) = -2.0/8.0;
	pixel(0, 2) = -1.0/8.0;

	// middle column
	pixel(1, 0) = 0.0;
	pixel(1, 1) = 0.0;
	pixel(1, 2) = 0.0;

	// right hand column
	pixel(2, 0) = 1.0/8.0;
	pixel(2, 1) = 2.0/8.0;
	pixel(2, 2) = 1.0/8.0;
}

const SobelOperatorX& SobelOperatorX::instance()
{
	return theInstance;
}

// SobelOperatorY stuff

SobelOperatorY SobelOperatorY::theInstance;

SobelOperatorY::SobelOperatorY()
:
FilterKernel3x3()
{
	// top row
	pixel(0, 0) = -1.0/8.0;
	pixel(1, 0) = -2.0/8.0;
	pixel(2, 0) = -1.0/8.0;

	// middle row
	pixel(0, 1) = 0.0;
	pixel(1, 1) = 0.0;
	pixel(2, 1) = 0.0;

	// right hand column
	pixel(0, 2) = 1.0/8.0;
	pixel(1, 2) = 2.0/8.0;
	pixel(2, 2) = 1.0/8.0;
}

const SobelOperatorY& SobelOperatorY::instance()
{
	return theInstance;
}

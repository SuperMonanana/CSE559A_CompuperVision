
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: pyramid.cpp                                                                          //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"




Pyramid::Pyramid(const Image& base, int levels)
:
Array<Image>(levels)
{
	//Image blurred;
	(*this)[0] = base;
	for (int i=1; i<levels; i++) {
		Image* current = &(*this)[i-1];
		Image* next = &(*this)[i];
		*next = Image((current->getWidth()+1)/2, (current->getHeight()+1)/2, base.getColors());
		current->resample(*next);
	}
}

Pyramid::Pyramid(int width, int height, int colors, int levels)
:
Array<Image>(levels)
{
	(*this)[0].resize(width, height, colors);
	for (int i=1; i<levels; i++) {
		Image* current = &(*this)[i-1];
		Image* next = &(*this)[i];
		next->resize((current->getWidth()+1)/2, (current->getHeight()+1)/2, current->getColors());
	}
}

void Pyramid::debugOutput(std::string filenamePattern) const
{
	for (int i=0; i<getSize(); i++) {
		std::string filename = Functions::filenameNumber(filenamePattern, i, getSize()-1);
		(*this)[i].saveTarga(filename);
	}
}



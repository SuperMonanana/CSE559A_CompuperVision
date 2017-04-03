
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: gaussian.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


// public

Gaussian::Gaussian(int size, double standard_deviation)
:
Vector(size)
{
	double mean = double(size-1)/2.0;
	for (int i=0; i<size; i++) {
		double entry = gaussian(double(i), mean, standard_deviation);
		Vector::operator[](i)=entry;
	}
}



double Gaussian::operator[](int i) const
{
	if (i>=getSize() || i<0) {
		return 0.0;
	}
	return Vector::operator[](i);
}

void Gaussian::output() const
{
	for (int i=0; i<(getSize()+1)/2; i++) {
		std::cout << (*this)[i] << std::endl;
	}
}



double Gaussian::gaussian(double x, double mean, double standard_deviation)
{
	double sqdev = standard_deviation*standard_deviation;
	double exponent = -((x-mean)*(x-mean)) / (2.0 * sqdev);
	return (1.0/sqrt(2.0*M_PI*sqdev)) * exp(exponent);
}






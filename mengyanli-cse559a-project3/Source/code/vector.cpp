
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: vector.cpp                                                                           //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"

Vector::Vector()
:
Array<double>()
{
	//empty
}

Vector::Vector(int size)
:
Array<double>(size)
{
	//empty
}

Vector& Vector::operator*=(double value)
{
	for (int i=0; i<getSize(); i++) {
		(*this)[i]*=value;
	}
	return *this;
}

Vector& Vector::operator/=(double value)
{
	for (int i=0; i<getSize(); i++) {
		(*this)[i]/=value;
	}
	return *this;
}

Vector& Vector::operator+=(double value)
{
	for (int i=0; i<getSize(); i++) {
		(*this)[i]+=value;
	}
	return *this;
}

Vector& Vector::operator+=(const Vector& rhs)
{
	if (rhs.getSize() != getSize()) {
		throw Error("+=: Arrays must be of same dimensions");
	}
	for (int i=0; i<getSize(); i++) {
		(*this)[i]+=rhs[i];
	}
	return *this;
}

void Vector::sub(const Vector& rhs, Vector& result) const
{
	if (rhs.getSize() != getSize()) {
		throw Error("Sub: Arrays must be of same dimensions");
	}
	result.resize(getSize());
	for (int i=0; i<getSize(); i++) {
		result[i]=(*this)[i]-rhs[i];
	}
}

void Vector::add(const Vector& rhs, Vector& result) const
{
	if (rhs.getSize() != getSize()) {
		throw Error("Add: Arrays must be of same dimensions");
	}
	result.resize(getSize());
	for (int i=0; i<getSize(); i++) {
		result[i]=(*this)[i]+rhs[i];
	}
}

double Vector::dot(const Vector& rhs) const
{
	if (rhs.getSize() != getSize()) {
		throw Error("Dot: Arrays must be of same dimensions");
	}
	double result=0.0;
	for (int i=0; i<getSize(); i++) {
		result+=((*this)[i]*rhs[i]);
	}
	return result;
}

double Vector::mse(const Vector& rhs) const
{
	if (rhs.getSize() != getSize()) {
		throw Error("Mse: Arrays must be of same dimensions");
	}
	double result=0.0;
	for (int i=0; i<getSize(); i++) {
		double error = ((*this)[i]-rhs[i]);
		result+=(error*error);
	}
	return result/getSize();
}

double Vector::mag() const
{
	double result=0.0;
	for (int i=0; i<getSize(); i++) {
		result+=((*this)[i]*(*this)[i]);
	}
	return sqrt(result);
}

double Vector::var() const
{
	double average=0.0;
	for (int i=0; i<getSize(); i++) {
		average+=(*this)[i];
	}
	average/=getSize();
	double squared_diff=0.0;
	for (int i=0; i<getSize(); i++) {
		double diff = (*this)[i]-average;
		squared_diff+=(diff*diff);
	}
	squared_diff/=getSize();
	return squared_diff;
}

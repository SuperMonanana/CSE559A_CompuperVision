

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: user.cpp                                                                             //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"

User::User()
:
Vector()
{
	//empty
}

User::User(std::string name, int size)
:
Vector(size),
name(name)
{
	//empty
}

void User::save(BinaryFileWriter& file) const
{
	file.write(name);
	Vector::save(file);
}

void User::load(BinaryFileReader& file)
{
	name=file.readString();
	Vector::load(file);
}

const std::string& User::getName() const
{
	return name;
}

void User::setName(const std::string& username)
{
	name=username;
}

bool User::operator<(const User& other) const
{
	return the_mse<other.the_mse;
}

void User::setMse(double value) const
{
	the_mse=value;
}

double User::getMse() const
{
	return the_mse;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: error.cpp                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


Error::Error(std::string message)
:
message(message)
{
	// empty
}

Error::Error(std::string message, std::string subs)
:
message(Functions::substitute(message, "%%", subs))
{
	// empty
}

Error::Error(std::string message, std::string subs1, std::string subs2)
:
message(Functions::substitute(Functions::substitute1(message, "%%", subs1), "%%", subs2))
{
	// empty
}

Error::Error(std::string message, std::string subs1, std::string subs2, std::string subs3)
:
message(Functions::substitute(Functions::substitute1(Functions::substitute1(message, "%%", subs1), "%%", subs2), "%%", subs3))
{
	// empty
}

Error::Error(std::string message, int subs)
:
message(Functions::substitute(message, "%%", Functions::toString(subs)))
{
	// empty
}

Error::Error(std::string message, int subs1, int subs2)
:
message(Functions::substitute(Functions::substitute1(message, "%%", Functions::toString(subs1)), "%%", Functions::toString(subs2)))
{
	// empty
}

Error::Error(std::string message, double subs)
:
message(Functions::substitute(message, "%%", Functions::toString(subs)))
{
	// empty
}

Error::Error(std::string message, double subs1, double subs2)
:
message(Functions::substitute(Functions::substitute1(message, "%%", Functions::toString(subs1)), "%%", Functions::toString(subs2)))
{
	// empty
}

const std::string& Error::getMessage() const
{
	return message;
}
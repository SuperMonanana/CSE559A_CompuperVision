
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: functions.cpp                                                                        //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


std::string Functions::substitute(std::string string, std::string marker, std::string subs)
{
	std::string::size_type pos=0;
	std::string::size_type length=marker.length();
	while(1) {
		pos = string.find(marker);
		if (pos<0 || pos>=string.length()) break;
		string.replace(pos, length, subs);
	}
	return string;
}

std::string Functions::substitute1(std::string string, std::string marker, std::string subs)
{
	std::string::size_type pos=0;
	std::string::size_type length=marker.length();
	pos = string.find(marker);
	if (pos>=0 && pos<string.length()) {
		string.replace(pos, length, subs);
	}
	return string;
}

std::string Functions::filenameNumber(std::string filenamePattern, int n, int max)
{
	std::string::size_type number_length = toString(max).length();
	return substitute(filenamePattern, "%%", toString(n, number_length));
}


std::string Functions::toString(int value)
{
	char buffer[32];
	sprintf(buffer, "%d", value);
	return std::string(buffer);
}

std::string Functions::toString(double value)
{
	char buffer[32];
	sprintf(buffer, "%f", value);
	return std::string(buffer);
}

std::string Functions::toString(int value, std::string::size_type length)
{
	char buffer[32];
	sprintf(buffer, "%d", value);
	std::string return_val(buffer);
	return_val = std::string(max(0, length-return_val.length()), '0') + return_val;
	return return_val;
}

void Functions::erase(std::string filename)
{
	std::cout << std::string(filename.length(), '\b');
}



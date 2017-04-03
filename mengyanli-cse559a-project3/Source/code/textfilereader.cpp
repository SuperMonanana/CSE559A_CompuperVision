
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: textfilereader.cpp                                                                   //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"


TextFileReader::TextFileReader(std::string filename)
:
FileReader(filename, std::ios_base::in)
{
	// empty
}

void TextFileReader::read(int& value)
{
	file >> value;
}

void TextFileReader::read(char& value)
{
	file >> value;
}

void TextFileReader::read(unsigned char& value)
{
	file >> value;
}

void TextFileReader::read(unsigned short& value)
{
	file >> value;
}

void TextFileReader::read(double& value)
{
	file >> value;
}

void TextFileReader::read(std::string& value)
{
	file >> value;
}
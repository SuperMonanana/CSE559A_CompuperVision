
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: binaryfilereader.cpp                                                                 //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


BinaryFileReader::BinaryFileReader(std::string filename)
:
FileReader(filename, std::ios_base::binary)
{
	// empty
}

void BinaryFileReader::read(int& value)
{
	file.read((char*)(&value), sizeof(value));
}

void BinaryFileReader::read(char& value)
{
	file >> value;
}

void BinaryFileReader::read(unsigned char& value)
{
	file.read((char*)(&value), sizeof(value));
}

void BinaryFileReader::read(unsigned short& value)
{
	file.read((char*)(&value), sizeof(value));
}

void BinaryFileReader::read(double& value)
{
	file.read((char*)(&value), sizeof(value));
}

void BinaryFileReader::read(std::string& value)
{
	value="";
	char c;
	while(true) {
		read(c);
		if (c=='\0' || file.eof()) break;
		value += c;
	}
}


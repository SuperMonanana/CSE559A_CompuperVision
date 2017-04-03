
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: filereader.cpp                                                                       //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


FileReader::FileReader(const std::string& filename, OPENMODE flags)
:
file(filename.c_str(), flags),
filename(filename)
{
	if (!file.is_open()) {
		throw Error("Can't open file '%%' for reading", filename);
	}
}

FileReader::~FileReader()
{
	file.close();
}

int FileReader::readInt()
{
	int result;
	read(result);
	return result;
}

char FileReader::readChar()
{
	char result;
	read(result);
	return result;

}

unsigned char FileReader::readUnsignedChar()
{
	unsigned char result;
	read(result);
	return result;

}

unsigned short FileReader::readUnsignedShort()
{
	unsigned short result;
	read(result);
	return result;

}

double FileReader::readDouble()
{
	double result;
	read(result);
	return result;
}

std::string FileReader::readString()
{
	std::string result;
	read(result);
	return result;
}

void FileReader::skip(int n)
{
	file.ignore(n);
}

const std::string& FileReader::getFilename()
{
	return filename;
}




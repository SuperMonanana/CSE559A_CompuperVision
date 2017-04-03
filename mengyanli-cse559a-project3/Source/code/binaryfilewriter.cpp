
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: binaryfilewriter.cpp                                                                 //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


BinaryFileWriter::BinaryFileWriter(std::string filename)
:
FileWriter(filename, std::ios_base::binary)
{
	// empty
}

void BinaryFileWriter::write(int value)
{
	file.write((char*)(&value), sizeof(value));
}

void BinaryFileWriter::write(char value)
{
	file << value;

}

void BinaryFileWriter::write(unsigned char value)
{
	file.write((char*)(&value), sizeof(value));

}

void BinaryFileWriter::write(unsigned short value)
{
	file.write((char*)(&value), sizeof(value));

}

void BinaryFileWriter::write(double value)
{
	file.write((char*)(&value), sizeof(value));

}

void BinaryFileWriter::write(std::string value)
{
	file.write(value.c_str(), int(value.length()));
	file << '\0';
}

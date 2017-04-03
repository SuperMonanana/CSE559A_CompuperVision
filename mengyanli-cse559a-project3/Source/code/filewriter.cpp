
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: filewriter.cpp                                                                       //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


FileWriter::FileWriter(const std::string& filename, OPENMODE flags)
:
file(filename.c_str(), flags),
filename(filename)
{
	if (!file.is_open()) {
		throw Error("Can't open file '%%' for writing", filename);
	}
}

FileWriter::~FileWriter()
{
	file.close();
}

const std::string& FileWriter::getFilename()
{
	return filename;
}


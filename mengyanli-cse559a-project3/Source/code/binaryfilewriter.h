
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: binaryfilewriter.h                                                                   //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class BinaryFileWriter : public FileWriter
{
public:
	BinaryFileWriter(std::string filename);
	virtual void write(int value);
	virtual void write(char value);
	virtual void write(unsigned char value);
	virtual void write(unsigned short value);
	virtual void write(double value);
	virtual void write(std::string value);
};

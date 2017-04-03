
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: binaryfilereader.h                                                                   //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class BinaryFileReader : public FileReader
{
public:
	BinaryFileReader(std::string filename);
	virtual void read(int& value);
	virtual void read(char& value);
	virtual void read(unsigned char& value);
	virtual void read(unsigned short& value);
	virtual void read(double& value);
	virtual void read(std::string& value);

};

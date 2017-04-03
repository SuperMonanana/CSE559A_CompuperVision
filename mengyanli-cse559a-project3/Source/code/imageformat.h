
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: imageformat.h                                                                        //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



class Image;

class ImageFormat
{
public:
	virtual void load(FileReader& file, Image& img)=0;
	virtual void load(std::string filename, Image& img);
	virtual void save(FileWriter& file, const Image& img)=0;
	virtual void save(std::string filename, const Image& img);
};
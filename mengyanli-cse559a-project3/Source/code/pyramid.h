
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: pyramid.h                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class Image;

class Pyramid : public Array<Image>
{
public:
	Pyramid(const Image& base, int levels);
	Pyramid(int width, int height, int colors, int levels);
	void debugOutput(std::string filenamePattern) const;
private:
	
};


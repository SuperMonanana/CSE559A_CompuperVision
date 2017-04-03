
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: functions.h                                                                          //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

namespace Functions {
	std::string substitute(std::string string, std::string marker, std::string subs);
	std::string substitute1(std::string string, std::string marker, std::string subs);
	std::string filenameNumber(std::string filenamePattern, int n, int max);
	std::string toString(int value);
	std::string toString(double value);
	std::string toString(int value, std::string::size_type length);
	void erase(std::string filename);

}
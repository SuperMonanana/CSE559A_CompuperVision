
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: error.h                                                                              //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



class Error
{
public:
	Error(std::string message);
	Error(std::string message, std::string subs);
	Error(std::string message, std::string subs1, std::string subs2);
	Error(std::string message, std::string subs1, std::string subs2, std::string subs3);
	Error(std::string message, int subs);
	Error(std::string message, int subs1, int subs2);
	Error(std::string message, double subs);
	Error(std::string message, double subs1, double subs2);
	const std::string& getMessage() const;

private:
	std::string message;
};




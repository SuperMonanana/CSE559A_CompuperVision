

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigenfaces.cpp                                                                       //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



// eigenfaces.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"




int main(int argc, char* argv[])
{
	try {
		Main main(argc, argv);
		main.run();
	}
	catch (Error e) {
		std::cout << std::endl << "Error: " << e.getMessage() << std::endl;
	}
	return 0;
}

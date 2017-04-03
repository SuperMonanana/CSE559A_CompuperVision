
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: commandlinearguments.cpp                                                             //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"


CommandLineArguments::CommandLineArguments(int argc, char** argv)
:
argc(argc),
argv(argv),
argIndex(1)
{
	// empty
}

std::string CommandLineArguments::getArg() const
{
	if (argIndex>=argc) {
		throw Error("Additional arguments expected");
	}
	return argv[argIndex];
}

std::string CommandLineArguments::nextArg()
{
	argIndex++;
	return getArg();
}

bool CommandLineArguments::moreArgs() const
{
	return ( argIndex < argc );
}




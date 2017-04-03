
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: commandlinearguments.h                                                               //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



class CommandLineArguments : public Arguments
{
public:
	CommandLineArguments(int argc, char** argv);
	virtual std::string getArg() const;
	virtual std::string nextArg();
	virtual bool moreArgs() const;

private:
	int argc;
	char** argv;
	int argIndex;
};
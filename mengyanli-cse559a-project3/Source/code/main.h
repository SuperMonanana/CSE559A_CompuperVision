

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: main.h                                                                               //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class Main
{
public:
	Main(int argc, char* argv[]);
	void run();
private:
	void usage() const;
	// command functions
	void computeEigenfaces();
	void projectFace();
	void constructUserbase();
	void isFace();
	void verifyFace();
	void recognizeFace();
	void findFace();
	void morphFaces();
	void invalidCommand();
	// end command functions
	CommandLineArguments args;
	std::string programName;

	// array of pointers to members
	void (Main::*commandFunctions[int(NUMBER_OF_COMMANDS)])();

};
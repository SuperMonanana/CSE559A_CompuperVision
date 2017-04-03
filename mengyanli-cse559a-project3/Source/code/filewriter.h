
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: filewriter.h                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef WIN32
#define OPENMODE std::ios_base::open_mode
#else
#define OPENMODE std::ios_base::openmode
#endif

class FileWriter
{
public:
	FileWriter(const std::string& filename, OPENMODE flags);
	~FileWriter();
	virtual void write(int value)=0;
	virtual void write(char value)=0;
	virtual void write(unsigned char value)=0;
	virtual void write(double value)=0;
	virtual void write(unsigned short value)=0;
	virtual void write(std::string value)=0;

	const std::string& getFilename();

protected:
	std::ofstream file;
	std::string filename;
private: // copying, assignment not allowed
	FileWriter(const FileWriter& writer) { /*empty*/ }
	FileWriter& operator=(const FileWriter& lhs) { /*empty*/ }
};


/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: filereader.h                                                                         //
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

class FileReader
{
public:
        FileReader(const std::string& filename, OPENMODE flags);
	~FileReader();
	virtual int readInt();
	virtual char readChar();
	virtual unsigned char readUnsignedChar();
	virtual unsigned short readUnsignedShort();
	virtual double readDouble();
	virtual std::string readString();
	virtual void read(int& value)=0;
	virtual void read(char& value)=0;
	virtual void read(unsigned char& value)=0;
	virtual void read(unsigned short& value)=0;
	virtual void read(double& value)=0;
	virtual void read(std::string& value)=0;


	void skip(int distance);
	const std::string& getFilename();
protected:
	std::ifstream file;
	std::string filename;
private: // copying, assignment not allowed
	FileReader(const FileReader& reader) { /*empty*/ }
	FileReader& operator=(const FileReader& lhs) { /*empty*/ }
};

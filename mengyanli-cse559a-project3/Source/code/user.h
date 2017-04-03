
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: user.h                                                                               //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

class User : public Vector
{
public:
	User();
	User(std::string name, int size);
	void save(BinaryFileWriter& file) const;
	void load(BinaryFileReader& file);
	const std::string& getName() const;
	void setName(const std::string& username);
	bool operator<(const User& other) const;
	void setMse(double value) const;
	double getMse() const;
private:
	std::string name;
	// just a temporary value for sorting purposes, so make it mutable
	mutable double the_mse;

};


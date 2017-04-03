
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: users.h                                                                              //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class Users : private std::list<User>
{
public:
	Users();
	void save(std::string filename) const;
	void save(BinaryFileWriter& file) const;
	void load(std::string filename);
	void load(BinaryFileReader& file);
	int getSize() const;
	void addUser(const User& user);
	const User& operator[](int index) const;
	const User& operator[](std::string) const;
	void sort();
private:
	std::map<std::string, const User*> map_string;
	std::map<int, const User*> map_int;

};
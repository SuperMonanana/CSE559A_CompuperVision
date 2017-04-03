
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: vector.h                                                                             //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

class Vector : public Array<double>
{
public:
	Vector();
	Vector(int size);

	Vector& operator*=(double value);
	Vector& operator/=(double value);
	Vector& operator+=(double value);
	Vector& operator+=(const Vector& rhs);
	void add(const Vector& rhs, Vector& result) const;
	void sub(const Vector& rhs, Vector& result) const;
	double dot(const Vector& rhs) const;
	double mse(const Vector& other) const;
	double mag() const;
	double var() const;
};

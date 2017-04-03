
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: filterkernel.h                                                                       //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class FilterKernel : private Image
{
protected:
	FilterKernel(int height, int width);
	double& pixel(int i, int j);
public:
	double pixel(int x, int y) const;
	int getHeight() const;
	int getWidth() const;
protected:
};



class FilterKernel3x3 : public FilterKernel
{
protected:
	FilterKernel3x3();
};

class FilterKernel5x5 : public FilterKernel
{
protected:
	FilterKernel5x5();
};

class BlurOperator3x3 : public FilterKernel3x3
{
public:
	static const BlurOperator3x3& instance();
private:
	BlurOperator3x3();
	static BlurOperator3x3 theInstance;
};

class BlurOperator5x5 : public FilterKernel5x5
{
public:
	static const BlurOperator5x5& instance();
private:
	BlurOperator5x5();
	static BlurOperator5x5 theInstance;
};

class SobelOperatorX : public FilterKernel3x3
{
public:
	static const SobelOperatorX& instance();
private:
	SobelOperatorX();
	static SobelOperatorX theInstance;
};

class SobelOperatorY : public FilterKernel3x3
{
public:
	static const SobelOperatorY& instance();
private:
	SobelOperatorY();
	static SobelOperatorY theInstance;
};

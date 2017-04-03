

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: face.h                                                                               //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class Face : public Image
{
public:
	Face(int width, int height);
	Face();
	void subimage(int x_low, int x_high, int y_low, int y_high, const Image& img, bool debug);
	void load(BinaryFileReader &file);
	void save(BinaryFileWriter &file) const;
	virtual void load(BinaryFileReader& file, ImageFormat& format);
private:

};


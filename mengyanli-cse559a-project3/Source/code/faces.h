

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: faces.h                                                                              //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

class EigFaces;

// for sorting eigenvectors
struct EigenVectorIndex
{
	double eigenvalue;
	int index;
	bool operator< (const EigenVectorIndex& rhs) const { return eigenvalue>rhs.eigenvalue; }
};

class Faces : public Array<Face>
{
public:
	Faces();
	Faces(int count, int width, int height);
	void load(BinaryFileReader& file);
	void load(std::string filename);
	void save(std::string filename) const;
	void save(BinaryFileWriter& file) const;
	void eigenFaces(EigFaces& results, int n) const;
	void output(std::string filepattern) const;
	int getWidth() const;
	int getHeight() const;
	void setWidth(int width);
	void setHeight(int height);
private:
	void sortEigenvalues(double *eigenvec, Array<int>& ordering) const;
protected:
	int width;
	int height;
	int vector_size;
	Face average_face;
};

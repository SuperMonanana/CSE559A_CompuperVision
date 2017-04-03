
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: image.h                                                                              //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


class ImageFormat;
class FilterKernel;

class Image : public Vector
{
public:
	Image();
	Image(int colors);
	Image(int height, int width);
	Image(int width, int height, int colors);
	inline double pixel(int i, int j, int c) const;
	inline double& pixel(int i, int j, int c);
	inline double sample(double i, double j, int c) const;
	int getWidth() const;
	int getHeight() const;
	int getColors() const;
	void resize(int new_width, int new_height);
	void resize(int new_width, int new_height, int new_colors);
	void crop(int x_min, int y_min, int x_max, int y_max, Image& result) const;
	void line(int x_min, int y_min, int x_max, int y_max, double r, double g, double b);
	void grayscale(Image& dest) const;
	void resample(Image& dest) const;
	double pixelFilter(const FilterKernel& kernel, int x, int y, int c) const;
	void imageFilter(Image& dest, const FilterKernel& kernel) const;

	void add(const Image& rhs, Image& result) const;
	void sub(const Image& rhs, Image& result) const;

	void normalize(double lo, double hi, Image& result) const;
	void saveTarga(std::string filename) const;
	void saveTarga(BinaryFileWriter& file) const;
	void save(std::string filename, ImageFormat& format) const;
	void save(BinaryFileWriter& file, ImageFormat& format) const;
	virtual void loadTarga(std::string filename);
	virtual void loadTarga(BinaryFileReader& file);
	virtual void load(std::string filename, ImageFormat& format);
	virtual void load(BinaryFileReader& file, ImageFormat& format);
protected:
	int width;
	int height;
	int colors;
};
double Image::pixel(int i, int j, int c) const
{
	if (i>=width || i<0 || j>=height || j<0 || c>=colors || c<0) {
		return(0.0);
	}

	return operator[]( colors * (j * width + i) + c);
}

double& Image::pixel(int i, int j, int c)
{
	static double dummy;
	if (i>=width || i<0 || j>=height || j<0 || c>=colors || c<0) {
		return(dummy=0.0);
	}

	return operator[]( colors * (j * width + i) + c);
}

double Image::sample(double i, double j, int color) const
{
	double a=pixel(int(i), int(j), color), b=pixel(int(i)+1, int(j), color), c=pixel(int(i), int(j)+1, color), d=pixel(int(i)+1, int(j)+1, color);
	double d1=i-int(i), d2=j-int(j);
	double d3=1.0-d1, d4=1.0-d2;
	double result= (d3*a+d1*b)*d4 + (d3*c+d1*d)*d2;
	return result;

}

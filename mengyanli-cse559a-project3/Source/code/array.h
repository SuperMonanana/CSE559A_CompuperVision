

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: array.h                                                                              //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


template <class T>
class Array : private std::vector<T>
{
public:
	Array();
	Array(int size);
	void save(BinaryFileWriter& file) const;
	void load(BinaryFileReader& file);
	inline int getSize() const;
	inline const T& operator[](int i) const;
	inline T& operator[](int i);
	void initialize();
	void initialize(const T& value);
	void resize(int new_size);
private:
	int size;

};

// public

template <class T>
Array<T>::Array()
:
std::vector<T>(),
size(0)
{
	// empty
}

template <class T>
Array<T>::Array(int size)
:
std::vector<T>(size),
size(size)
{
	// emtpy
}

template <class T>
void Array<T>::save(BinaryFileWriter& file) const
{
	file.write(size);
	for (int i=0; i<size; i++) {
		file.write((*this)[i]);
	}
}

template <class T>
void Array<T>::load(BinaryFileReader& file)
{
	resize(file.readInt());
	for (int i=0; i<size; i++) {
		T value;
		file.read(value);
		(*this)[i]=value;
	}
}

template <class T>
int Array<T>::getSize() const
{
	return size;
}

template <class T>
const T& Array<T>::operator[](int i) const
{
	if (i<0 || i>=size) {
		throw Error("Array index out of range");
	}
	return std::vector<T>::operator [](i);
}

template <class T>
T& Array<T>::operator[](int i)
{
	if (i<0 || i>=size) {
		throw Error("Array index out of range");
	}
	return std::vector<T>::operator [](i);
}

template <class T>
void Array<T>::initialize()
{
	initialize(0);
}

template <class T>
void Array<T>::initialize(const T& value)
{
	for (int i=0; i<size; i++) {
		(*this)[i] = value;
	}
}

template <class T>
void Array<T>::resize(int new_size)
{
	if (new_size != size) {
		std::vector<T>::resize(new_size);
		size=new_size;
	}
}




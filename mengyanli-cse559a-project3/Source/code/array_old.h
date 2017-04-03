

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
class Array
{
public:
	Array();
	Array(int size);
	~Array();
	Array(const Array<T>& array);
	Array<T>& operator=(const Array<T>& rhs);
	void save(BinaryFileWriter& file) const;
	void load(BinaryFileReader& file);
	inline int getSize() const;
	inline const T& operator[](int i) const;
	inline T& operator[](int i);
	Array<T>& operator*=(const T& value);
	Array<T>& operator/=(const T& value);
	Array<T>& operator+=(const T& value);
	Array<T>& operator+=(const Array<T>& rhs);
	//Array<T> operator-(const Array<T>& rhs) const;
	//Array<T> operator+(const Array<T>& rhs) const;
	void add(const Array<T>& rhs, Array<T>& result) const;
	void sub(const Array<T>& rhs, Array<T>& result) const;
	T dot(const Array<T>& rhs) const;
	T mse(const Array<T>& other) const;
	void initialize();
	void initialize(const T& value);
	void resize(int new_size);
protected:
	void reallocate(int new_size);
private:
	void copy(const Array<T>& array);
	void clean();
	void allocate(int new_size);
	int size;
	T* data;
};

// public

template <class T>
Array<T>::Array()
:
data(NULL),
size(0)
{
	// empty
}

template <class T>
Array<T>::Array(int size)
:
data(NULL)
{
	allocate(size);
}

template <class T>
Array<T>::~Array()
{
	clean();
}

template <class T>
Array<T>::Array(const Array<T>& array)
:
data(NULL)
{
	copy(array);
}

template <class T>
Array<T>& Array<T>::operator=(const Array<T>& rhs)
{
	if (&rhs!=this) {
		copy(rhs);
	}
	return *this;
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
	return data[i];
}

template <class T>
T& Array<T>::operator[](int i)
{
	if (i<0 || i>=size) {
		throw Error("Array index out of range");
	}
	return data[i];
}

template <class T>
Array<T>& Array<T>::operator*=(const T& value)
{
	for (int i=0; i<size; i++) {
		data[i]*=value;
	}
	return *this;
}

template <class T>
Array<T>& Array<T>::operator/=(const T& value)
{
	for (int i=0; i<size; i++) {
		data[i]/=value;
	}
	return *this;
}

template <class T>
Array<T>& Array<T>::operator+=(const T& value)
{
	for (int i=0; i<size; i++) {
		data[i]+=value;
	}
	return *this;
}

template <class T>
Array<T>& Array<T>::operator+=(const Array<T>& rhs)
{
	if (rhs.size != size) {
		throw Error("+=: Arrays must be of same dimensions");
	}
	for (int i=0; i<size; i++) {
		data[i]+=rhs.data[i];
	}
	return *this;
}

template <class T>
void Array<T>::sub(const Array<T>& rhs, Array<T>& result) const
{
	if (rhs.size != size) {
		throw Error("Sub: Arrays must be of same dimensions");
	}
	result.resize(size);
	for (int i=0; i<size; i++) {
		result[i]=(*this)[i]-rhs[i];
	}
}

template <class T>
void Array<T>::add(const Array<T>& rhs, Array<T>& result) const
{
	if (rhs.size != size) {
		throw Error("Add: Arrays must be of same dimensions");
	}
	result.resize(size);
	for (int i=0; i<size; i++) {
		result[i]=(*this)[i]+rhs[i];
	}
}

template <class T>
T Array<T>::dot(const Array<T>& rhs) const
{
	if (rhs.size != size) {
		throw Error("Dot: Arrays must be of same dimensions");
	}
	T result=0.0;
	for (int i=0; i<size; i++) {
		result+=((*this)[i]*rhs[i]);
	}
	return result;
}

template <class T>
T Array<T>::mse(const Array<T>& rhs) const
{
	if (rhs.size != size) {
		throw Error("Mse: Arrays must be of same dimensions");
	}
	T result=0.0;
	for (int i=0; i<size; i++) {
		T error = ((*this)[i]-rhs[i]);
		result+=(error*error);
	}
	return result/size;
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
		data[i] = value;
	}
}

template <class T>
void Array<T>::resize(int new_size)
{
	allocate(new_size);
}

// protected

template <class T>
void Array<T>::reallocate(int new_size)
{
	T* new_data = new T[new_size];
	for (int i=0; i<size && i<new_size; i++) {
		new_data[i] = data[i];
	}
	clean();
	size=new_size;
	data=new_data;
}

// private

template <class T>
void Array<T>::copy(const Array<T>& array)
{
	allocate(array.size);
	for (int i=0; i<size; i++) {
		data[i] = array.data[i];
	}
}

template <class T>
void Array<T>::clean()
{
	if (data) {
		delete [] data;
	}
	size = 0;
}

template <class T>
void Array<T>::allocate(int new_size)
{
	if (!data || size!=new_size) {
		clean();
		size=new_size;
		data = new T[size];
	}
}


typedef Array<double> Vector;

# Makefile for ImageLib

IMAGELIB=libImage.a
IMAGELIB_OBJS=Convert.o Convolve.o FileIO.o Image.o ImageProc.o Pyramid.o \
		RefCntMem.o Transform.o WarpImage.o

CC=g++
CPPFLAGS=-Wall -O3

all: $(IMAGELIB)

$(IMAGELIB): 
	g++ -c -o Convert.o Convert.cpp
	g++ -c -o Convolve.o Convolve.cpp
	g++ -c -o FileIO.o FileIO.cpp
	g++ -c -o Image.o Image.cpp
	g++ -c -o ImageProc.o ImageProc.cpp
	g++ -c -o Pyramid.o Pyramid.cpp
	g++ -c -o RefCntMem.o RefCntMem.cpp
	g++ -c -o Transform.o Transform.cpp
	g++ -c -o WarpImage.o WarpImage.cpp
	ar r $(IMAGELIB) $(IMAGELIB_OBJS)

clean:
	rm -f *.o *~ $(IMAGELIB)

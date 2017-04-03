

/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: stdafx.h                                                                             //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once


#ifdef WIN32
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#define OPENMODE std::ios_base::open_mode
#else
#define OPENMODE std::ios_base::openmode
#endif



#include <stdio.h>
//#include <tchar.h>

// reference additional headers your program requires here

#define _USE_MATH_DEFINES

enum Direction {
	LEFT,
	RIGHT
};



#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "functions.h"

#include "error.h"
#include "filereader.h"
#include "binaryfilereader.h"
#include "textfilereader.h"
#include "filewriter.h"
#include "binaryfilewriter.h"
#include "array.h"
#include "vector.h"
#include "image.h"
#include "face.h"
#include "user.h"
#include "users.h"
#include "faces.h"
#include "eigfaces.h"
#include "filterkernel.h"
#include "imageformat.h"
#include "targa.h"
#include "pyramid.h"
#include "gaussian.h"
//#include "commandmap.h"
#include "arguments.h"
#include "commandlinearguments.h"
#include "main.h"


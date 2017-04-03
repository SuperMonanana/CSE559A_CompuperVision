
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: arguments.cpp                                                                        //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include "stdafx.h"



CommandMap Arguments::commandMap;

Command Arguments::nextCommand()
{
	return toCommand(nextArg());
}

Command Arguments::getCommand() const
{
	return toCommand(getArg());
}

int Arguments::nextPositiveInt()
{
	int value = atoi(nextArg().c_str());
	if (value < 0) throw Error("Expected non-negative integer, got %%", value);
	return value;
}

int Arguments::getPositiveInt()
{
	int value = atoi(getArg().c_str());
	if (value < 0) throw Error("Expected non-negative integer, got %%", value);
	return value;
}

double Arguments::nextFloat()
{
	double value = atof(nextArg().c_str());
	return value;
}

double Arguments::getFloat()
{
	double value = atof(getArg().c_str());
	return value;
}

double Arguments::nextPositiveFloat()
{
	double value = atof(nextArg().c_str());
	if (value < 0.0) throw Error("Expected non-negative float, got %%", value);
	return value;
}

double Arguments::getPositiveFloat()
{
	double value = atof(getArg().c_str());
	if (value < 0.0) throw Error("Expected non-negative float, got %%", value);
	return value;
}

std::string Arguments::nextExpect(std::string expected)
{
	std::string value = nextArg();
	if (value!=expected) {
		throw Error("Expected '%%', got '%%'", expected, value);
	}
	return value;
}

std::string Arguments::getExpect(std::string expected)
{
	std::string value = getArg();
	if (value!=expected) {
		throw Error("Expected '%%', got '%%'", expected, value);
	}
	return value;
}

std::string Arguments::nextExpect(std::string expected1, std::string expected2)
{
	std::string value = nextArg();
	if (value!=expected1 && value!=expected2) {
		throw Error("Expected '%%' or '%%', got '%%'", expected1, expected2, value);
	}
	return value;
}

std::string Arguments::getExpect(std::string expected1, std::string expected2)
{
	std::string value = getArg();
	if (value!=expected1 && value!=expected2) {
		throw Error("Expected '%%' or '%%', got '%%'", expected1, expected2, value);
	}
	return value;
}

Command Arguments::toCommand(std::string string)
{
	return commandMap.getCommand(string);
}

CommandMap::CommandMap()
{
	map["--eigenfaces"]		= COMPUTE_EIGENFACES;
	map["--projectface"]    = PROJECT_FACE;
	map["--constructuserbase"] = CONSTRUCT_USERBASE;
	map["--isface"]         = IS_FACE;
	map["--verifyface"]     = VERIFY_FACE;
	map["--recognizeface"]  = RECOGNIZE_FACE;
	map["--findface"]       = FIND_FACE;
	map["--morphfaces"]     = MORPH_FACES;

}

Command CommandMap::getCommand(const std::string& name)
{
	if (map.find(name)!=map.end()) {
		return map[name];
	}
	else return INVALID_COMMAND;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: arguments.h                                                                          //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



enum Command
{
	INVALID_COMMAND,
	COMPUTE_EIGENFACES,
	PROJECT_FACE,
	CONSTRUCT_USERBASE,
	IS_FACE,
	VERIFY_FACE,
	RECOGNIZE_FACE,
	FIND_FACE,
	MORPH_FACES,
	NUMBER_OF_COMMANDS
};

class CommandMap
{
public:
	CommandMap();
	Command getCommand(const std::string& name);

private:
	std::map<std::string, Command> map;
};


class Arguments
{
public:
	virtual Command nextCommand();
	virtual Command getCommand() const;
	virtual int nextPositiveInt();
	virtual int getPositiveInt();
	virtual double nextFloat();
	virtual double getFloat();
	virtual double nextPositiveFloat();
	virtual double getPositiveFloat();
	virtual std::string nextExpect(std::string expected);
	virtual std::string getExpect(std::string expected);
	virtual std::string nextExpect(std::string expected1, std::string expected2);
	virtual std::string getExpect(std::string expected1, std::string expected2);
	virtual std::string nextArg()=0;
	virtual std::string getArg() const=0;
	virtual bool moreArgs() const=0;
protected:
	static Command toCommand(std::string string);

	static CommandMap commandMap;
};

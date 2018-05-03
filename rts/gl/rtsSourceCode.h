#ifndef RTSSOURCECODE_H
#define RTSSOURCECODE_H

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

///This class defines generic source code that can be loaded from text files.  It is primarily used by the rts_glShaderProgram class for GLSL programming.

class rtsSourceCode
{
public:
	vector<string> source;			//the actual source code
	void clear()					///<Clears any current source code from the class.
	{
		source.clear();
	}
	void LoadSource(const char* filename)	///<Loads source code from a specified file.
	{
		ifstream infile;		//create an input file
		infile.open(filename);	//load the specified file
	
		if(!infile.is_open())	//if the file is not open, exit
		{
			return;
		}
		source.clear();			//remove any previous code

		while(!infile.eof())
		{
			string current_line;		
			getline(infile, current_line);
			current_line += '\n';
			source.push_back(current_line);
		}
	}
	rtsSourceCode(const char* filename)	///<Constructor creates the class and loads source code from the specified file.
	{
		LoadSource(filename);
	}
	rtsSourceCode(){}						///<Constructor creates a blank class.
	rtsSourceCode& operator+=(const rtsSourceCode& rhs)
	{
		int lines = rhs.source.size();
		for(int l=0; l<lines; l++)
			source.push_back(rhs.source[l]);
		return *this;
	}
	rtsSourceCode& operator+=(const string& rhs)
	{
		source.push_back(rhs);
		return *this;
	}
	void ConsoleOut()						///<Sends the source code to the standard output.
	{
		unsigned int lines = source.size();
		for(unsigned int l = 0; l<lines; l++)
			cout<<l<<":  "<<source[l];
	}
};

#endif

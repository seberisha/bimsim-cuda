#ifndef RTS_GLSHADERS
#define RTS_GLSHADERS

#include <GL/glew.h>
//#include "windows.h"
#include <GL/gl.h>
#include "rtsSourceCode.h"

class rts_glShaderObject
{
private:
	void init()
	{
		id = 0;
		compiled = false;
		type = GL_FRAGMENT_SHADER;
	}
public:
	bool compiled;
	GLenum type;
	rtsSourceCode source;
	GLuint id;
	string log;

	rts_glShaderObject(GLenum type, const char* filename)
	{
		init();					//initialize the shader
		SetType(type);	//set the shader type
		LoadSource(filename);	//load the source code
	}
	rts_glShaderObject(GLenum type, rtsSourceCode sourceCode)
	{
		init();					//initialize the shader
		SetType(type);	//set the shader type
		source = sourceCode;
	}
	rts_glShaderObject()
	{
		init();
	}
	rts_glShaderObject(GLenum type)
	{
		init();
		SetType(type);
	}
	void LoadSource(const char* filename)
	{
		source = rtsSourceCode(filename);	//get the shader source code

	}
	void SetType(GLenum type)
	{
		if(id != 0)					//if a shader currently exists, delete it
		{
			glDeleteShader(id);
			id = 0;
		}
		type = type;
		id = glCreateShader(type);		//create a shader object
		if(id == 0)						//if a shader was not created, log an error
		{
			log = "Error getting shader ID from OpenGL";
			return;
		}
	}
	void UploadSource()
	{
		//create the structure for the shader source code
		GLsizei count = source.source.size();
		GLchar** code_string = new GLchar*[count];
		GLint* length = new GLint[count];
		for(int l = 0; l<count; l++)	//for each line of code
		{
			length[l] = source.source[l].size();
			code_string[l] = new GLchar[length[l]];	//copy the string into a new structure
			source.source[l].copy(code_string[l], (unsigned int)length[l]);

		}
		glShaderSource(id, count, (const GLchar**)code_string, length);		//attach the shader source
	}
	void Compile()
	{
		/*
		This function compiles the shader source code, records any errors to a log, and sets the compiled flag.
		*/
		//send the source code to the GPU
		UploadSource();

		//compile the shader
		glCompileShader(id);												//compile the shader
		GLint compile_status;
		glGetShaderiv(id, GL_COMPILE_STATUS, &compile_status);				//get the compile status
		if(compile_status != GL_TRUE)	//if there was an error
		{
			GLchar buffer[1000];		//create a log buffer
			GLsizei length;
			glGetShaderInfoLog(id, 1000, &length, buffer);	//get the log
			log = buffer;
			compiled = false;
		}
		else
			compiled = true;

	}
	void PrintLog()
	{
		cout<<log;
		if(log.size() != 0) cout<<endl;
	}
	void Clean(){if(id != 0) glDeleteShader(id);}
};



#endif

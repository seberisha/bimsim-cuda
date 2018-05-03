#ifndef RTS_GLSHADERUNIFORM_H
#define RTS_GLSHADERUNIFORM_H

#include "CHECK_OPENGL_ERROR.h"
#include <GL/glew.h>
#include <string>

using namespace std;

enum rtsUniformEnum {RTS_FLOAT, RTS_INT, RTS_BOOL, RTS_FLOAT_MATRIX};

///This class stores a single uniform variable for GLSL and is designed to be used by the rts_glShaderProgram class.
struct rts_glShaderUniform
{
public:
	string name;		//the name of the variable
	GLint location;		//the location in the program
	void* p_value;		//pointer to the global data representing the value in main memory
	GLenum type;		//variable type (float, int, vec2, etc.)
	//rtsUniformEnum rts_type;	//type of variable in rts format
	//unsigned int num;		//the number of values required by the variable (1 for float, 2 for vec2, etc.)
	string log;

	//void convert_type(GLenum gl_type);		//converts the OpenGL data type to something useful for rts
	void submit_to_gpu()
	{
		if(location < 0)
			return;
		if(p_value == NULL)
		{
			cout<<"Error in uniform address: "<<name<<endl;
			return;
		}
	

		CHECK_OPENGL_ERROR
		switch(type)
		{
		case GL_FLOAT:
			glUniform1fv(location, 1, (float*)p_value);
			break;
		case GL_FLOAT_VEC2:
			glUniform2fv(location, 1, (float*)p_value);
			break;
		case GL_FLOAT_VEC3:
			glUniform3fv(location, 1, (float*)p_value);
			break;
		case GL_FLOAT_VEC4:
			glUniform4fv(location, 1, (float*)p_value);
			break;
		case GL_INT:
			glUniform1iv(location, 1, (int*)p_value);
			break;
		case GL_INT_VEC2:
			glUniform2iv(location, 1, (int*)p_value);
			break;
		case GL_INT_VEC3:
			glUniform3iv(location, 1, (int*)p_value);
			break;
		case GL_INT_VEC4:
			glUniform4iv(location, 1, (int*)p_value);
			break;
		case GL_BOOL:
			glUniform1iv(location, 1, (int*)p_value);
			break;
		case GL_BOOL_VEC2:
			glUniform2iv(location, 1, (int*)p_value);
			break;
		case GL_BOOL_VEC3:
			glUniform3iv(location, 1, (int*)p_value);
			break;
		case GL_BOOL_VEC4:
			glUniform4iv(location, 1, (int*)p_value);
			break;
		case GL_FLOAT_MAT2:
			glUniformMatrix2fv(location, 1, GL_FALSE, (float*)p_value);
			break;
		case GL_FLOAT_MAT3:
			glUniformMatrix3fv(location, 1, GL_FALSE, (float*)p_value);
			break;
		case GL_FLOAT_MAT4:
			glUniformMatrix4fv(location, 1, GL_FALSE, (float*)p_value);
			break;
		case GL_SAMPLER_1D:
		case GL_SAMPLER_2D:
		case GL_SAMPLER_3D:
		case GL_SAMPLER_CUBE:
		case GL_SAMPLER_1D_SHADOW:
		case GL_SAMPLER_2D_SHADOW:
		default:
			glUniform1iv(location, 1, (int*)p_value);
			break;
		}
		CHECK_OPENGL_ERROR
	}
	rts_glShaderUniform()
	{
		location = -1;
		p_value = NULL;
	}
};



#endif

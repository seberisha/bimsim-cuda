#ifndef RTS_GLTEXTUREMAP_H
#define RTS_GLTEXTUREMAP_H

//#include <GL/glew.h>
#include "rts/math/vector.h"
#include "rts/gl/error.h"
#include <stdlib.h>

namespace rts{

///This class stores an OpenGL texture map and is used by rts_glShaderProgram.
class glTexture
{
private:
	void get_type()			//guesses the texture type based on the size
	{
		if(size[1] == 0)
			texture_type = GL_TEXTURE_1D;
		else if(size[2] == 0)
			texture_type = GL_TEXTURE_2D;
		else
			texture_type = GL_TEXTURE_3D;
	}
	void set_wrapping()		//set the texture wrapping based on the dimensions
	{
		CHECK_OPENGL_ERROR
		switch(texture_type)
		{
		case GL_TEXTURE_3D:
			glTexParameteri(texture_type, GL_TEXTURE_WRAP_R_EXT, GL_REPEAT);
		case GL_TEXTURE_2D:
			glTexParameteri(texture_type, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
		case GL_TEXTURE_1D:
			glTexParameteri(texture_type, GL_TEXTURE_WRAP_S, GL_REPEAT);
			break;
		case GL_TEXTURE_RECTANGLE_ARB:
			glTexParameteri(texture_type, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(texture_type, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			break;

		default:
			break;
		}
		CHECK_OPENGL_ERROR
	}
	//void set_bits(GLvoid* bits);
public:
	vector<GLsizei, 3> size;		//vector representing the size of the texture
	GLuint name;				//texture name assigned by OpenGL
	GLenum texture_type;				//1D, 2D, 3D
	GLint internal_format;	//number of components (ex. 4 for RGBA)
	GLenum pixel_format;		//type of data (RGBA, LUMINANCE)
	GLenum data_type;			//data type of the bits (float, int, etc.)

	//constructor
	glTexture()
	{
		name = 0;
	}
	glTexture(GLvoid *bits,
			   GLenum type = GL_TEXTURE_2D,
			   GLsizei width = 256,
			   GLsizei height = 256,
			   GLsizei depth = 0,
			   GLint internalformat = 1,
			   GLenum format = GL_LUMINANCE,
			   GLenum datatype = GL_UNSIGNED_BYTE,
			   GLint interpolation = GL_LINEAR)
    {
        init(bits, type, width, height, depth, internalformat, format, datatype, interpolation);
    }

	void begin()
	{
		glEnable(texture_type);
		CHECK_OPENGL_ERROR
		glBindTexture(texture_type, name);
		CHECK_OPENGL_ERROR
	}
	void end()
	{
		glDisable(texture_type);
		CHECK_OPENGL_ERROR
	}

	///Creates an OpenGL texture map. This function requires basic information about the texture map as well as a pointer to the bit data describing the texture.
	void init(GLvoid *bits,
			   GLenum type = GL_TEXTURE_2D,
			   GLsizei width = 256,
			   GLsizei height = 256,
			   GLsizei depth = 0,
			   GLint internalformat = 1,
			   GLenum format = GL_LUMINANCE,
			   GLenum datatype = GL_UNSIGNED_BYTE,
			   GLint interpolation = GL_LINEAR)
	{
		CHECK_OPENGL_ERROR
		if(name != 0)
			glDeleteTextures(1, &name);


		CHECK_OPENGL_ERROR
		if(datatype == GL_FLOAT)
		{
			glPixelStorei(GL_PACK_ALIGNMENT, 4);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 4);				//I honestly don't know what this does but it fixes problems
		}
		else if(datatype == GL_UNSIGNED_BYTE)
		{
			//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			//glPixelStorei(GL_PACK_ALIGNMENT, 1);
		}
		else if(datatype == GL_UNSIGNED_SHORT)
		{
			//glPixelStorei(GL_UNPACK_ALIGNMENT, 2);
			//glPixelStorei(GL_PACK_ALIGNMENT, 2);
		}
		CHECK_OPENGL_ERROR
		glGenTextures(1, &name);							//get the texture name from OpenGL
		//cout<<"OpenGL Name: "<<name<<endl;
		CHECK_OPENGL_ERROR
		size = vector<GLsizei, 3>(width, height, depth);		//assign the texture size
		//get_type();											//guess the type based on the size
		texture_type = type;						//set the type of texture
		glEnable(texture_type);						//enable the texture map
		CHECK_OPENGL_ERROR
		glBindTexture(texture_type, name);							//bind the texture for editing
		CHECK_OPENGL_ERROR
		set_wrapping();										//set the texture wrapping parameters
		CHECK_OPENGL_ERROR
		glTexParameteri(texture_type, GL_TEXTURE_MAG_FILTER, interpolation);		//set filtering
		CHECK_OPENGL_ERROR
		glTexParameteri(texture_type, GL_TEXTURE_MIN_FILTER, interpolation);
		CHECK_OPENGL_ERROR
		internal_format = internalformat;					//set the number of components per pixel
		pixel_format = format;								//set the pixel format
		data_type = datatype;									//set the data type
		SetBits(bits);										//send the bits to the OpenGL driver
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);	//replace the specified vertex color
		CHECK_OPENGL_ERROR
		glDisable(texture_type);

		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	}
	void Clean()
	{
		if(name != 0)
		{
			glDeleteTextures(1, &name);
			CHECK_OPENGL_ERROR
			name = 0;
		}
	}
	void SetBits(GLvoid *bits)
	{
		glEnable(texture_type);						//enable the texture map
		CHECK_OPENGL_ERROR
		glBindTexture(texture_type, name);
		CHECK_OPENGL_ERROR

		switch(texture_type)
		{
		case GL_TEXTURE_3D:
			glTexImage3D(texture_type, 0, internal_format, size[0], size[1], size[2], 0, pixel_format, data_type, bits);
			CHECK_OPENGL_ERROR
			break;
		case GL_TEXTURE_2D:
		case GL_TEXTURE_RECTANGLE_ARB:
			glTexImage2D(texture_type, 0, internal_format, size[0], size[1], 0, pixel_format, data_type, bits);
			CHECK_OPENGL_ERROR
			break;
		case GL_TEXTURE_1D:
			glTexImage1D(texture_type, 0, internal_format, size[0], 0, pixel_format, data_type, bits);
			CHECK_OPENGL_ERROR
			break;
		default:
			//glTexImage2D(texture_type, 0, internal_format, size.x, size.y, 0, pixel_format, data_type, bits);
			break;
		}
		CHECK_OPENGL_ERROR
	}
	void ResetBits(GLvoid *bits)
	{
		glEnable(texture_type);						//enable the texture map
		CHECK_OPENGL_ERROR
		glBindTexture(texture_type, name);
		CHECK_OPENGL_ERROR

		switch(texture_type)
		{
		case GL_TEXTURE_3D:
			//glTexImage3D(texture_type, 0, internal_format, size.x, size.y, size.z, 0, pixel_format, data_type, bits);
			break;
		case GL_TEXTURE_2D:
		case GL_TEXTURE_RECTANGLE_ARB:
			glTexSubImage2D(texture_type, 0, 0, 0, size[0], size[1], pixel_format, data_type, bits);
			CHECK_OPENGL_ERROR
			break;
		case GL_TEXTURE_1D:
			//glTexImage1D(texture_type, 0, internal_format, size.x, 0, pixel_format, data_type, bits);
			break;
		default:
			//glTexImage2D(texture_type, 0, internal_format, size.x, size.y, 0, pixel_format, data_type, bits);
			break;
		}
		glDisable(texture_type);
		CHECK_OPENGL_ERROR
	}
	void* GetBits(GLenum format, GLenum type)
	{
		//returns the texture data

		int components;
		switch(format)
		{
		case GL_RED:
		case GL_GREEN:
		case GL_BLUE:
		case GL_ALPHA:
		case GL_LUMINANCE:
			components = 1;
			break;
		case GL_LUMINANCE_ALPHA:
			components = 2;
			break;
		case GL_RGB:
		case GL_BGR:
			components = 3;
			break;
		case GL_RGBA:
		case GL_BGRA:
			components = 4;
			break;
		}

		int type_size;
		switch(type)
		{
		case GL_UNSIGNED_BYTE:
		case GL_BYTE:
			type_size = sizeof(char);
			break;
		case GL_UNSIGNED_SHORT:
		case GL_SHORT:
			type_size = sizeof(short);
			break;
		case GL_UNSIGNED_INT:
		case GL_INT:
			type_size = sizeof(int);
			break;
		case GL_FLOAT:
			type_size = sizeof(float);
			break;
		}

		//allocate memory for the texture
		void* result = malloc(components*type_size * size[0] * size[1]);

		begin();
		glGetTexImage(texture_type, 0, format, type, result);

		CHECK_OPENGL_ERROR
		end();


		return result;

	}
};

}

#define RTS_UNKNOWN				0

#endif

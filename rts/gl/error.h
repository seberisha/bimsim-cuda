#ifndef RTS_OPENGL_ERROR
#define RTS_OPENGL_ERROR

#include <stdio.h>
#include <GL/gl.h>
#include <GL/glu.h>

#define CHECK_OPENGL_ERROR \
{ GLenum error; \
   while ( (error = glGetError()) != GL_NO_ERROR) { \
   printf( "OpenGL ERROR: %s\nCHECK POINT: %s (line %d)\n", gluErrorString(error), __FILE__, __LINE__ ); \
   } \
}

#endif

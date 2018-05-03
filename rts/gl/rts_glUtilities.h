#ifndef RTS_GLUTILITIES_H
#define RTS_GLUTILITIES_H


#define CHECK_OPENGL_ERROR \
{ GLenum error; \
   while ( (error = glGetError()) != GL_NO_ERROR) { \
   printf( "OpenGL ERROR: %s\nCHECK POINT: %s (line %d)\n", gluErrorString(error), __FILE__, __LINE__ ); \
   } \
} 

#endif

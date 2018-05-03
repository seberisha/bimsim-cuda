#ifndef DATATYPES_H
#define DATATYPES_H

#include <math.h>

#define PRECISION_SINGLE

#ifdef PRECISION_SINGLE
typedef float ptype;
#elif defined PRECISION_DOUBLE
typedef double ptype;
#endif

#define BLOCK	256
#define SQRT_BLOCK 8

#define PI	3.14159f

//a very small number
#define EPSILON		0.00001f

//CUDA hybrid code - complex class should run on both the CPU and GPU


typedef ptype fieldPoint;

extern bool verbose;

//hybrid GPU/CPU complex data typ
#include "rts/math/complex.h"
#include "rts/math/vector.h"
#include "rts/math/point.h"
#include "rts/math/quad.h"

typedef rts::complex<ptype> bsComplex;
typedef rts::vector<ptype, 3> bsVector;
typedef rts::point<ptype, 3> bsPoint;
typedef rts::quad<ptype, 3> bsRect;


#endif

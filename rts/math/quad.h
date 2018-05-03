#ifndef RTS_RECT_H
#define RTS_RECT_H

//enable CUDA_CALLABLE macro
#include "cuda_callable.h"
#include "rts/math/vector.h"
#include "rts/math/point.h"
#include "rts/math/triangle.h"
#include "rts/math/quaternion.h"
#include <iostream>
#include <algorithm>

namespace rts{

//template for a quadangle class in ND space
template <class T, int N>
struct quad
{
	/*
		C------------------>O
		^                   ^
		|                   |
		Y                   |
		|                   |
		|                   |
		A---------X-------->B
	*/

	/*T A[N];
	T B[N];
	T C[N];*/

        rts::point<T, N>A;
	rts::vector<T, N> X;
	rts::vector<T, N> Y;


	CUDA_CALLABLE quad()
	{

	}

	CUDA_CALLABLE quad(point<T, N> a, point<T, N> b, point<T, N> c)
	{

		A = a;
		X = b - a;
		Y = c - a;

	}

    /****************************************************************
    Constructor - create a quad from two points and a normal
    ****************************************************************/
	CUDA_CALLABLE quad(rts::point<T, N> pMin, rts::point<T, N> pMax, rts::vector<T, N> normal)
	{

        //assign the corner point
        A = pMin;

        //compute the vector from pMin to pMax
        rts::vector<T, 3> v0;
        v0 = pMax - pMin;
       // std::cout<<"v0==="<<v0.toStr()<<std::endl;


        //compute the cross product of A and the plane normal
        rts::vector<T, 3> v1;
        v1 = v0.cross(normal);
      // std::cout<<"normal==="<<normal.toStr()<<std::endl;
       // std::cout<<"v1==="<<v1.toStr()<<std::endl;

        //calculate point B
        rts::point<T, 3> B;
        B = A + v0 * 0.5 + v1 * 0.5;
     //   std::cout<<"B==="<<B.toStr()<<std::endl;
        //calculate rtsPoint C
        rts::point<T, 3> C;
        C = A  + v0 * 0.5 - v1 * 0.5;
   //     std::cout<<"C==="<<C.toStr()<<std::endl;
        //calculate X and Y
        X = B - A;
        Y = C - A;
	}

	/*******************************************************************
	Constructor - create a quad from a position, normal, and rotation
	*******************************************************************/
	CUDA_CALLABLE quad(rts::point<T, N> c, rts::vector<T, N> normal, T width, T height, T theta)
	{

        //compute the X direction - start along world-space X
        Y = rts::vector<T, N>(0, 1, 0);
        if(Y == normal)
            Y = rts::vector<T, N>(0, 0, 1);

        X = Y.cross(normal).norm();

     //   std::cout<<X<<std::endl;

        //rotate the X axis by theta radians
        rts::quaternion<T> q;
        q.CreateRotation(theta, normal);
        X = q.toMatrix3() * X;
        Y = normal.cross(X);

        //normalize everything
        X = X.norm();
        Y = Y.norm();

        //scale to match the quad width and height
        X = X * width;
        Y = Y * height;

        //set the corner of the plane
        A = c - X * 0.5 - Y * 0.5;

     //   std::cout<<X<<std::endl;
	}

	/*******************************************
	Return the normal for the quad
	*******************************************/
	CUDA_CALLABLE rts::vector<T, N> n()
	{
        return (X.cross(Y)).norm();
	}

        CUDA_CALLABLE rts::point<T, N> p(T a, T b)
	{
		rts::point<T, N> result;
		//given the two parameters a, b = [0 1], returns the position in world space
//                std::cout<<"A==="<<A<<std::endl;
//                std::cout<<"X==="<<X<<std::endl;
//                std::cout<<"a==="<<a<<std::endl;
//                std::cout<<"Y==="<<Y<<std::endl;
//                std::cout<<"b==="<<b<<std::endl;

		result = A + X * a + Y * b;

		return result;
	}

	CUDA_CALLABLE rts::point<T, N> operator()(T a, T b)
	{
		return p(a, b);
	}

	std::string toStr()
	{
		std::stringstream ss;

		ss<<"A = "<<A<<std::endl;
		ss<<"B = "<<A + X<<std::endl;
		ss<<"C = "<<A + X + Y<<std::endl;
		ss<<"D = "<<A + Y<<std::endl;

        return ss.str();

	}

	CUDA_CALLABLE quad<T, N> operator*(T rhs)
	{
		//scales the plane by a scalar value

		//compute the center point
		rts::point<T, N> c = A + X*0.5 + Y*0.5;

		//create the new quadangle
		quad<T, N> result;
		result.X = X * rhs;
		result.Y = Y * rhs;
		result.A = c - result.X*0.5 - result.Y*0.5;

		return result;

	}

	CUDA_CALLABLE T dist(point<T, N> p)
	{
        //compute the distance between a point and this quad

        //first break the quad up into two triangles
        triangle<T, N> T0(A, A+X, A+Y);
        triangle<T, N> T1(A+X+Y, A+X, A+Y);


        ptype d0 = T0.dist(p);
        ptype d1 = T1.dist(p);

        if(d0 < d1)
            return d0;
        else
            return d1;
	}

	CUDA_CALLABLE T dist_max(point<T, N> p)
	{
        T da = (A - p).len();
        T db = (A+X - p).len();
        T dc = (A+Y - p).len();
        T dd = (A+X+Y - p).len();

        return std::max( da, std::max(db, std::max(dc, dd) ) );
	}
};

}	//end namespace rts

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, rts::quad<T, N> R)
{
    os<<R.toStr();
    return os;
}


#endif

#ifndef RTS_TRIANGLE_H
#define RTS_TRIANGLE_H

//enable CUDA_CALLABLE macro
#include "cuda_callable.h"
#include "rts/math/vector.h"
#include "rts/math/point.h"
#include <iostream>

namespace rts{

template <class T, int N>
struct triangle
{
    /*
        A------>B
        |      /
        |     /
        |    /
        |   /
        |  /
        | /
        C
    */
    private:

    point<T, N> A;
    point<T, N> B;
    point<T, N> C;

    CUDA_CALLABLE point<T, N> _p(T s, T t)
    {
        //This function returns the point specified by p = A + s(B-A) + t(C-A)
        vector<T, N> E0 = B-A;
        vector<T, N> E1 = C-A;

        return A + s*E0 + t*E1;
    }


    public:



    CUDA_CALLABLE triangle()
	{

	}

	CUDA_CALLABLE triangle(point<T, N> a, point<T, N> b, point<T, N> c)
	{
		A = a;
		B = b;
		C = c;
	}

	CUDA_CALLABLE rts::point<T, N> operator()(T s, T t)
	{
        return _p(s, t);
	}

	CUDA_CALLABLE point<T, N> nearest(point<T, N> p)
	{
        //comptue the distance between a point and this triangle
        //  This code is adapted from: http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

        vector<T, N> E0 = B-A;
        vector<T, N> E1 = C-A;
        vector<T, N> D = A - p;

        T a = E0.dot(E0);
        T b = E0.dot(E1);
        T c = E1.dot(E1);
        T d = E0.dot(D);
        T e = E1.dot(D);
        //T f = D.dot(D);

        T det = a*c - b*b;
        T s = b*e - c*d;
        T t = b*d - a*e;

        /*std::cout<<"E0: "<<E0<<std::endl;
        std::cout<<"E1: "<<E1<<std::endl;
        std::cout<<"a: "<<a<<std::endl;
        std::cout<<"b: "<<b<<std::endl;
        std::cout<<"c: "<<c<<std::endl;
        std::cout<<"d: "<<d<<std::endl;
        std::cout<<"e: "<<e<<std::endl;
        std::cout<<"f: "<<f<<std::endl;
        std::cout<<"det: "<<det<<std::endl;
        std::cout<<"s: "<<s<<std::endl;
        std::cout<<"t: "<<t<<std::endl;*/


        if( s+t <= det)
        {
            if(s < 0)
            {
                if(t < 0)
                {
                    //region 4
                    //std::cout<<"Region 4"<<std::endl;
                    s = 0;
                    t = 0;
                    //done?
                }
                else
                {
                    //region 3
                    //std::cout<<"Region 3"<<std::endl;
                    s=0;
                    t = ( e >= 0 ? 0 : ( -e >= c ? 1 : -e/c ) );
                    //done
                }
            }
            else if(t < 0)
            {
                //region 5
                //std::cout<<"Region 5"<<std::endl;
                s = ( d >= 0 ? 0 : ( -d >= a ? 1 : -d/a ) );
                t = 0;
                //done
            }
            else
            {
                //region 0
                //std::cout<<"Region 0"<<std::endl;
                T invDet = (ptype)1.0/det;
                s *= invDet;
                t *= invDet;
                //done
            }
        }
        else
        {
            if(s < 0)
            {
                //region 2
                //std::cout<<"Region 2"<<std::endl;
                s = 0;
                t = 1;
                //done?

            }
            else if(t < 0)
            {
                //region 6
                //std::cout<<"Region 6"<<std::endl;
                s = 1;
                t = 0;
                //done?
            }
            else
            {
                //region 1
                //std::cout<<"Region 1"<<std::endl;
                T numer = c + e - b - d;
                if( numer <= 0 )
                    s = 0;
                else
                {
                    T denom = a - 2 * b + c;
                    s = ( numer >= denom ? 1 : numer/denom );
                }
                t = 1 - s;
                //done
            }
        }

        //std::cout<<"s: "<<s<<std::endl;
        //std::cout<<"t: "<<t<<std::endl;

        //std::cout<<"p: "<<_p(s, t)<<std::endl;

		return _p(s, t);

	}

	CUDA_CALLABLE T dist(point<T, N> p)
	{
        point<T, N> n = nearest(p);

        return (p - n).len();
	}
};

}

#endif

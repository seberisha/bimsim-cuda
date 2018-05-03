#ifndef RTS_VECTOR_H
#define RTS_VECTOR_H

#include <iostream>
#include <cmath>
#include <sstream>
//#include "rts/point.h"
#include "cuda_callable.h"

#include <stdio.h>

namespace rts
{



template <class T, int N>
struct vector
{
	T v[N];

	CUDA_CALLABLE vector()
	{
		//memset(v, 0, sizeof(T) * N);
		for(int i=0; i<N; i++)
			v[i] = 0;
	}

	//efficiency constructor, makes construction easier for 1D-4D vectors
	CUDA_CALLABLE vector(T x, T y = (T)0.0, T z = (T)0.0, T w = (T)0.0)
	{
		if(N >= 1)
			v[0] = x;
		if(N >= 2)
			v[1] = y;
		if(N >= 3)
			v[2] = z;
		if(N >= 4)
			v[3] = w;
	}

	CUDA_CALLABLE vector(const T(&data)[N])
	{
		memcpy(v, data, sizeof(T) * N);
	}

	CUDA_CALLABLE T len()
	{
        //compute and return the vector length
        T sum_sq = (T)0;
        for(int i=0; i<N; i++)
        {
            //printf("\t\t v[i]: %f\n", v[i] );
            sum_sq += v[i] * v[i];
        }
        return std::sqrt(sum_sq);

	}

	CUDA_CALLABLE vector<T, N> cart2sph()
	{
		//convert the vector from cartesian to spherical coordinates
		//x, y, z -> r, theta, phi (where theta = 0 to 2*pi)

                vector<T, N                                     > sph;
		sph[0] = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
                sph[1] = std::atan2(v[1], v[0]);
                sph[2] = std::acos(v[2] / sph[0]);

                return sph;
	}

	CUDA_CALLABLE vector<T, N> sph2cart()
	{
		//convert the vector from cartesian to spherical coordinates
		//r, theta, phi -> x, y, z (where theta = 0 to 2*pi)

		vector<T, N> cart;
		cart[0] = v[0] * std::cos(v[1]) * std::sin(v[2]);
		cart[1] = v[0] * std::sin(v[1]) * std::sin(v[2]);
		cart[2] = v[0] * std::cos(v[2]);

		return cart;
	}

	CUDA_CALLABLE vector<T, N> norm()
	{
        //compute and return the vector norm
        vector<T, N> result;

        //compute the vector length
        T l = len();

        //normalize
        for(int i=0; i<N; i++)
        {
            result.v[i] = v[i] / l;
        }

        return result;
	}

	CUDA_CALLABLE vector<T, 3> cross(vector<T, 3> rhs)
	{
		vector<T, 3> result;

		//compute the cross product (only valid for 3D vectors)
		result[0] = v[1] * rhs[2] - v[2] * rhs[1];
		result[1] = v[2] * rhs[0] - v[0] * rhs[2];
		result[2] = v[0] * rhs[1] - v[1] * rhs[0];

		return result;
	}

    CUDA_CALLABLE T dot(vector<T, N> rhs)
    {
        T result = (T)0;

        for(int i=0; i<N; i++)
            result += v[i] * rhs.v[i];

        return result;

    }

	//arithmetic
	CUDA_CALLABLE vector<T, N> operator+(vector<T, N> rhs)
	{
        vector<T, N> result;

        for(int i=0; i<N; i++)
            result.v[i] = v[i] + rhs.v[i];

        return result;
	}
	CUDA_CALLABLE vector<T, N> operator-(vector<T, N> rhs)
	{
        vector<T, N> result;

        for(int i=0; i<N; i++)
            result.v[i] = v[i] - rhs.v[i];

        return result;
	}
	CUDA_CALLABLE vector<T, N> operator*(T rhs)
	{
        vector<T, N> result;

        for(int i=0; i<N; i++)
            result.v[i] = v[i] * rhs;

        return result;
	}
	CUDA_CALLABLE vector<T, N> operator/(T rhs)
	{
        vector<T, N> result;

        for(int i=0; i<N; i++)
            result.v[i] = v[i] / rhs;

        return result;
	}

	CUDA_CALLABLE bool operator==(vector<T, N> rhs)
	{
        if ( (rhs.v[0] == v[0]) && (rhs.v[1] == v[1]) && (rhs.v[2] == v[2]) )
            return true;

        return false;
	}

	std::string toStr()
	{
		std::stringstream ss;

		ss<<"[";
		for(int i=0; i<N; i++)
		{
			ss<<v[i];
			if(i != N-1)
				ss<<", ";
		}
		ss<<"]";

		return ss.str();
	}

	//bracket operator
	CUDA_CALLABLE T& operator[](int i)
	{
        return v[i];
    }

};


}	//end namespace rts

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, rts::vector<T, N> v)
{
    os<<v.toStr();
    return os;
}

//arithmetic operators
template <typename T, int N>
CUDA_CALLABLE rts::vector<T, N> operator-(rts::vector<T, N> v)
{
    rts::vector<T, N> r;

    //negate the vector
    for(int i=0; i<N; i++)
        r.v[i] = -v.v[i];

    return r;
}

template <typename T, int N>
CUDA_CALLABLE rts::vector<T, N> operator*(T lhs, rts::vector<T, N> rhs)
{
    rts::vector<T, N> r;

    return rhs * lhs;
}

//#if __GNUC__ > 3 && __GNUC_MINOR__ > 7
//template<class T, int N> using rtsVector = rts::vector<T, N>;
//#endif

#endif

#ifndef RTS_rtsPoint_H
#define RTS_rtsPoint_H

#include "rts/math/vector.h"
#include <string.h>
#include "cuda_callable.h"

namespace rts
{

template <class T, int N>
struct point
{
	T p[N];

	CUDA_CALLABLE point()
	{

	}

	//efficiency constructor, makes construction easier for 1D-4D vectors
	CUDA_CALLABLE point(T x, T y = (T)0.0, T z = (T)0.0, T w = (T)0.0)
	{
		if(N >= 1)
			p[0] = x;
		if(N >= 2)
			p[1] = y;
		if(N >= 3)
			p[2] = z;
		//if(N >= 4)
		//	p[3] = w;
	}

	//arithmetic operators
	CUDA_CALLABLE rts::point<T, N> operator+(vector<T, N> v)
	{
        rts::point<T, N> r;

        //calculate the position of the resulting point
        for(int i=0; i<N; i++)
            r.p[i] = p[i] + v.v[i];

        return r;
	}
	CUDA_CALLABLE rts::point<T, N> operator-(vector<T, N> v)
	{
        rts::point<T, N> r;

        //calculate the position of the resulting point
        for(int i=0; i<N; i++)
            r.p[i] = p[i] - v.v[i];

        return r;
	}
	CUDA_CALLABLE vector<T, N> operator-(point<T, N> rhs)
	{
        vector<T, N> r;

        //calculate the position of the resulting point
        for(int i=0; i<N; i++)
            r.v[i] = p[i] - rhs.p[i];

        return r;
	}
	CUDA_CALLABLE rts::point<T, N> operator*(T rhs)
	{
        rts::point<T, N> r;

        //calculate the position of the resulting point
        for(int i=0; i<N; i++)
            r.p[i] = p[i] * rhs;

        return r;
	}

	CUDA_CALLABLE point(const T(&data)[N])
	{
		memcpy(p, data, sizeof(T) * N);
	}

	std::string toStr()
	{
		std::stringstream ss;

		ss<<"(";
		for(int i=0; i<N; i++)
		{
			ss<<p[i];
			if(i != N-1)
				ss<<", ";
		}
		ss<<")";

		return ss.str();
	}

	//bracket operator
	CUDA_CALLABLE T& operator[](int i)
	{
        return p[i];
    }

};

}	//end namespace rts

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, rts::point<T, N> p)
{
    os<<p.toStr();
    return os;
}

//arithmetic
template <typename T, int N>
CUDA_CALLABLE rts::point<T, N> operator*(T lhs, rts::point<T, N> rhs)
{
    rts::point<T, N> r;

    return rhs * lhs;
}

//#if __GNUC__ > 3 && __GNUC_MINOR__ > 7
//template<class T, int N> using rtsPoint = rts::point<T, N>;
//#endif

#endif

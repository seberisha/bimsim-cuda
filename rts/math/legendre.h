#ifndef RTS_LEGENDRE_H
#define RTS_LEGENDRE_H

#include "cuda_callable.h"

namespace rts{

template <typename T>
CUDA_CALLABLE void init_legendre(T x, T& P0, T& P1)
{
	//compute the first two Legendre polynomials
	P0 = 1;
	P1 = x;
}

template <typename T>
CUDA_CALLABLE void shift_legendre(int n, T x, T& P0, T& P1)
{
	//compute the next (order n) Legendre polynomial
	T Pnew = ( (2 * n - 1) * x * P1 - (n-1) * P0 ) / n;

	//shift and add the new value to the array
	P0 = P1;
	P1 = Pnew;
}

template <typename T>
CUDA_CALLABLE void legendre(int n, T x, T* P)
{
    P[0] = 1;

    if(n >= 1)
        P[1] = x;

    for(int i=2; i<=n; i++)
    {
        P[i] = ( (2 * i - 1) * x * P[i-1] - (i-1) * P[i-2] ) / i;
    }

}

}


#endif

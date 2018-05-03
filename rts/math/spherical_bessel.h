#ifndef RTS_SBESSEL_H
#define RTS_SBESSEL_H
#include <math.h>


namespace rts{

#define RTS_BESSEL_CONVERGENCE_MIN		0.0145
#define RTS_BESSEL_CONVERGENCE_MAX		0.4
#define RTS_BESSEL_MAXIMUM_FLOAT		-1e33

template <typename T>
CUDA_CALLABLE void sbesselj(int n, complex<T> x, complex<T>* j)
{
    //compute the first bessel function
    if(n >= 0)
        j[0] = sin(x) / x;

    //compute the second bessel function
    if(n >= 1)
        j[1] = j[0] / x - cos(x) / x;

    //use the recurrence relation to compute the rest
    for(int i = 2; i <= n; i++)
    {
        j[i] = ( (2 * i - 1) / x ) * j[i-1] - j[i-2];
    }

    //if x = 0, deal with the degenerate case
    /*if( isnan(j[0].r) )
    {
        j[0] = (T)1.0;
        for(int i = 1; i<=n; i++)
            j[i] = (T)0.0;
    }*/
}

template <typename T>
CUDA_CALLABLE void sbessely(int n, complex<T> x, complex<T>* y)
{
    //compute the first bessel function
    if(n >= 0)
        y[0] = -cos(x) / x;

    //compute the second bessel function
    if(n >= 1)
        y[1] = y[0] / x - sin(x) / x;

    //use the recurrence relation to compute the rest
    for(int i = 2; i <= n; i++)
    {
        y[i] = ( (2 * i - 1) / x ) * y[i-1] - y[i-2];
    }

}

//spherical Hankel functions of the first kind
template <typename T>
CUDA_CALLABLE void sbesselh1(int n, complex<T> x, complex<T>* h)
{
    //compute j_0 and j_1
    complex<T> j[2];
    sbesselj(1, x, j);

    //compute y_0 and y_1
    complex<T> y[2];
    sbessely(1, x, y);

    //compute the first-order Hhankel function
    if(n >= 0)
        h[0] = j[0] + y[0].imul();

    //compute the second bessel function
    if(n >= 1)
        h[1] = j[1] + y[1].imul();

    //use the recurrence relation to compute the rest
    for(int i = 2; i <= n; i++)
    {
        h[i] = ( (2 * i - 1) / x ) * h[i-1] - h[i-2];
    }
}

template <typename T>
CUDA_CALLABLE void init_sbesselj(T x, T* j)
{
	//compute the first 2 bessel functions
	j[0] = sin(x) / x;

	j[1] = j[0] / x - cos(x) / x;
}

template <typename T>
CUDA_CALLABLE void init_sbessely(T x, T* y)
{
	//compute the first 2 bessel functions
	y[0] = -cos(x) / x;

	y[1] = y[0] / x - sin(x) / x;
}

template <typename T>
CUDA_CALLABLE void shift_sbesselj(int n, T x, T* b)//, T stability = 1.4)
{

	T bnew;

	//compute the next (order n) Bessel function
	bnew = ((2 * n - 1) * b[1])/x - b[0];

	//if(n > stability*x)
	if(n > real(x))
		if(real(bnew) < RTS_BESSEL_CONVERGENCE_MIN || real(bnew) > RTS_BESSEL_CONVERGENCE_MAX)
			bnew = 0.0;

	//shift and add the new value to the array
	b[0] = b[1];
	b[1] = bnew;
}

template <typename T>
CUDA_CALLABLE void shift_sbessely(int n, T x, T* b)//, T stability = 1.4)
{

	T bnew;

	//compute the next (order n) Bessel function
	bnew = ((2 * n - 1) * b[1])/x - b[0];

	if(bnew < RTS_BESSEL_MAXIMUM_FLOAT ||
	   (n > x && bnew > 0))
	{
		bnew = (T)0;
		b[1] = (T)0;
	}


	//shift and add the new value to the array
	b[0] = b[1];
	b[1] = bnew;
}



}   //end namespace rts



#endif

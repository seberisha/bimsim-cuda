/*RTS Complex number class.  This class is CUDA compatible,
and can therefore be used in CUDA code and on CUDA devices.
*/


#ifndef RTS_COMPLEX
#define RTS_COMPLEX

#include "cuda_callable.h"
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

namespace rts
{

template <class T>
struct complex
{
    T r, i;

    //default constructor
    CUDA_CALLABLE complex()
    {
        r = 0.0;
		i = 0.0;
    }

	//access methods
	CUDA_CALLABLE T real()
	{
		return r;
	}

	CUDA_CALLABLE T real(T r_val)
	{
		r = r_val;
		return r_val;
	}

	CUDA_CALLABLE T imag()
	{
		return i;
	}
	CUDA_CALLABLE T imag(T i_val)
	{
		i = i_val;
		return i_val;
	}

    //constructor when given real and imaginary values
    CUDA_CALLABLE complex(T r, T i)
    {
        this->r = r;
        this->i = i;
    }

    //return the current value multiplied by i
    CUDA_CALLABLE complex<T> imul()
    {
        complex<T> result;
        result.r = -i;
        result.i = r;

        return result;
    }

	//ARITHMETIC OPERATORS--------------------

    //binary + operator (returns the result of adding two complex values)
    CUDA_CALLABLE complex<T> operator+ (const complex<T> rhs)
    {
        complex<T> result;
        result.r = r + rhs.r;
        result.i = i + rhs.i;
        return result;
    }

	CUDA_CALLABLE complex<T> operator+ (const T rhs)
    {
        complex<T> result;
        result.r = r + rhs;
        result.i = i;
        return result;
    }

    //binary - operator (returns the result of adding two complex values)
    CUDA_CALLABLE complex<T> operator- (const complex<T> rhs)
    {
        complex<T> result;
        result.r = r - rhs.r;
        result.i = i - rhs.i;
        return result;
    }

    //binary - operator (returns the result of adding two complex values)
    CUDA_CALLABLE complex<T> operator- (const T rhs)
    {
        complex<T> result;
        result.r = r - rhs;
        result.i = i;
        return result;
    }

    //binary MULTIPLICATION operators (returns the result of multiplying complex values)
    CUDA_CALLABLE complex<T> operator* (const complex<T> rhs)
    {
        complex<T> result;
        result.r = r * rhs.r - i * rhs.i;
        result.i = r * rhs.i + i * rhs.r;
        return result;
    }
    CUDA_CALLABLE complex<T> operator* (const T rhs)
    {
        return complex<T>(r * rhs, i * rhs);
    }

    //binary DIVISION operators (returns the result of dividing complex values)
    CUDA_CALLABLE complex<T> operator/ (const complex<T> rhs)
    {
        complex<T> result;
        T denom = rhs.r * rhs.r + rhs.i * rhs.i;
        result.r = (r * rhs.r + i * rhs.i) / denom;
        result.i = (- r * rhs.i + i * rhs.r) / denom;

        return result;
    }
    CUDA_CALLABLE complex<T> operator/ (const T rhs)
    {
        return complex<T>(r / rhs, i / rhs);
    }

    //ASSIGNMENT operators-----------------------------------
    CUDA_CALLABLE complex<T> & operator=(const complex<T> &rhs)
    {
        //check for self-assignment
        if(this != &rhs)
        {
            this->r = rhs.r;
            this->i = rhs.i;
        }
        return *this;
    }
    CUDA_CALLABLE complex<T> & operator=(const T &rhs)
    {
        this->r = rhs;
        this->i = 0;

		return *this;
    }

    //arithmetic assignment operators
    CUDA_CALLABLE complex<T> operator+=(const complex<T> &rhs)
    {
		*this = *this + rhs;
        return *this;
    }
    CUDA_CALLABLE complex<T> operator+=(const T &rhs)
    {
		*this = *this + rhs;
        return *this;
    }

    CUDA_CALLABLE complex<T> operator*=(const complex<T> &rhs)
    {
		*this = *this * rhs;
        return *this;
    }
	CUDA_CALLABLE complex<T> operator*=(const T &rhs)
    {
		*this = *this * rhs;
        return *this;
    }
	//divide and assign
	CUDA_CALLABLE complex<T> operator/=(const complex<T> &rhs)
    {
		*this = *this / rhs;
        return *this;
    }
    CUDA_CALLABLE complex<T> operator/=(const T &rhs)
    {
		*this = *this / rhs;
        return *this;
    }

    //absolute value operator (returns the absolute value of the complex number)
	CUDA_CALLABLE T abs()
	{
		return std::sqrt(r * r + i * i);
	}

	CUDA_CALLABLE complex<T> log()
	{
        complex<T> result;
        result.r = std::log(std::sqrt(r * r + i * i));
        result.i = std::atan2(i, r);


        return result;
	}

	CUDA_CALLABLE complex<T> exp()
	{
        complex<T> result;

        T e_r = std::exp(r);
        result.r = e_r * std::cos(i);
        result.i = e_r * std::sin(i);

        return result;
	}

	/*CUDA_CALLABLE complex<T> pow(int y)
	{

        return pow((double)y);
	}*/

	CUDA_CALLABLE complex<T> pow(T y)
	{
        complex<T> result;

        result = log() * y;

        return result.exp();
	}

	CUDA_CALLABLE complex<T> sqrt()
	{
		complex<T> result;

		//convert to polar coordinates
		T a = std::sqrt(r*r + i*i);
		T theta = std::atan2(i, r);

		//find the square root
		T a_p = std::sqrt(a);
		T theta_p = theta/2.0;

		//convert back to cartesian coordinates
		result.r = a_p * std::cos(theta_p);
		result.i = a_p * std::sin(theta_p);

		return result;
	}

	std::string toStr()
	{
		std::stringstream ss;
		ss<<"("<<r<<","<<i<<")";

		return ss.str();
	}

	//COMPARISON operators
	CUDA_CALLABLE bool operator==(complex<T> rhs)
	{
        if(r == rhs.r && i == rhs.i)
            return true;
        return false;
    }

    CUDA_CALLABLE bool operator==(T rhs)
	{
        if(r == rhs && i == (T)0.0)
            return true;
        return false;
    }

};

}	//end RTS namespace

//addition
template<typename T>
CUDA_CALLABLE static rts::complex<T> operator+(const double a, const rts::complex<T> b)
{
    return rts::complex<T>(a + b.r, b.i);
}

//subtraction with a real value
template<typename T>
CUDA_CALLABLE static rts::complex<T> operator-(const double a, const rts::complex<T> b)
{
    return rts::complex<T>(a - b.r, -b.i);
}

//minus sign
template<typename T>
CUDA_CALLABLE static rts::complex<T> operator-(const rts::complex<T> &rhs)
{
    return rts::complex<T>(-rhs.r, -rhs.i);
}

//multiply a T value by a complex value
template<typename T>
CUDA_CALLABLE static rts::complex<T> operator*(const double a, const rts::complex<T> b)
{
    return rts::complex<T>((T)a * b.r, (T)a * b.i);
}

//divide a T value by a complex value
template<typename T>
CUDA_CALLABLE static rts::complex<T> operator/(const double a, const rts::complex<T> b)
{
    //return complex<T>(a * b.r, a * b.i);
    rts::complex<T> result;

    T denom = b.r * b.r + b.i * b.i;

    result.r = (a * b.r) / denom;
    result.i = -(a * b.i) / denom;

    return result;
}

//POW function
/*template<typename T>
CUDA_CALLABLE static complex<T> pow(complex<T> x, int y)
{
	return x.pow(y);
}*/

template<typename T>
CUDA_CALLABLE static rts::complex<T> pow(rts::complex<T> x, T y)
{
	return x.pow(y);
}

//log function
template<typename T>
CUDA_CALLABLE static rts::complex<T> log(rts::complex<T> x)
{
	return x.log();
}

//exp function
template<typename T>
CUDA_CALLABLE static rts::complex<T> exp(rts::complex<T> x)
{
	return x.exp();
}

//sqrt function
template<typename T>
CUDA_CALLABLE static rts::complex<T> sqrt(rts::complex<T> x)
{
	return x.sqrt();
}


template <typename T>
CUDA_CALLABLE static T abs(rts::complex<T> a)
{
    return a.abs();
}

template <typename T>
CUDA_CALLABLE static T real(rts::complex<T> a)
{
    return a.r;
}

//template <typename T>
CUDA_CALLABLE static float real(float a)
{
    return a;
}

template <typename T>
CUDA_CALLABLE static T imag(rts::complex<T> a)
{
    return a.i;
}

//trigonometric functions
template<class A>
CUDA_CALLABLE rts::complex<A> sin(const rts::complex<A> x)
{
	rts::complex<A> result;
	result.r = std::sin(x.r) * std::cosh(x.i);
	result.i = std::cos(x.r) * std::sinh(x.i);

	return result;
}

template<class A>
CUDA_CALLABLE rts::complex<A> cos(const rts::complex<A> x)
{
	rts::complex<A> result;
	result.r = std::cos(x.r) * std::cosh(x.i);
	result.i = -(std::sin(x.r) * std::sinh(x.i));

	return result;
}


template<class A>
std::ostream& operator<<(std::ostream& os, rts::complex<A> x)
{
    os<<x.toStr();
    return os;
}

//#if __GNUC__ > 3 && __GNUC_MINOR__ > 7
//template<class T> using rtsComplex = rts::complex<T>;
//#endif



#endif

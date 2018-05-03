#include "rts/math/matrix.h"

#ifndef RTS_QUATERNION_H
#define RTS_QUATERNION_H

namespace rts{

template<typename T>
class quaternion
{
public:
	T w;
	T x;
	T y;
	T z;

	void normalize();
	void CreateRotation(T theta, T axis_x, T axis_y, T axis_z);
	void CreateRotation(T theta, vector<T, 3> axis);
	quaternion<T> operator*(quaternion<T> &rhs);
	matrix<T, 3> toMatrix3();
	matrix<T, 4> toMatrix4();


	quaternion();
	quaternion(T w, T x, T y, T z);

};

template<typename T>
void quaternion<T>::normalize()
{
	double length=sqrt(w*w + x*x + y*y + z*z);
	w=w/length;
	x=x/length;
	y=y/length;
	z=z/length;
}

template<typename T>
void quaternion<T>::CreateRotation(T theta, T axis_x, T axis_y, T axis_z)
{
	//assign the given Euler rotation to this quaternion
	w = (T)cos(theta/2.0);
	x = axis_x*(T)sin(theta/2.0);
	y = axis_y*(T)sin(theta/2.0);
	z = axis_z*(T)sin(theta/2.0);
}

template<typename T>
void quaternion<T>::CreateRotation(T theta, vector<T, 3> axis)
{
	CreateRotation(theta, axis[0], axis[1], axis[2]);
}

template<typename T>
quaternion<T> quaternion<T>::operator *(quaternion<T> &param)
{
	float A, B, C, D, E, F, G, H;


	A = (w + x)*(param.w + param.x);
	B = (z - y)*(param.y - param.z);
	C = (w - x)*(param.y + param.z);
	D = (y + z)*(param.w - param.x);
	E = (x + z)*(param.x + param.y);
	F = (x - z)*(param.x - param.y);
	G = (w + y)*(param.w - param.z);
	H = (w - y)*(param.w + param.z);

	quaternion<T> result;
	result.w = B + (-E - F + G + H) /2;
	result.x = A - (E + F + G + H)/2;
	result.y = C + (E - F + G - H)/2;
	result.z = D + (E - F - G + H)/2;

	return result;
}

template<typename T>
matrix<T, 3> quaternion<T>::toMatrix3()
{
	matrix<T, 3> result;


    T wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;


    // calculate coefficients
    x2 = x + x; y2 = y + y;
    z2 = z + z;
    xx = x * x2; xy = x * y2; xz = x * z2;
    yy = y * y2; yz = y * z2; zz = z * z2;
    wx = w * x2; wy = w * y2; wz = w * z2;

	result(0, 0) = (T)1.0 - (yy + zz);
	result(0, 1) = xy - wz;

	result(0, 2) = xz + wy;

	result(1, 0) = xy + wz;
	result(1, 1) = (T)1.0 - (xx + zz);

	result(1, 2) = yz - wx;

	result(2, 0) = xz - wy;
	result(2, 1) = yz + wx;

	result(2, 2) = (T)1.0 - (xx + yy);

	return result;
}

template<typename T>
matrix<T, 4> quaternion<T>::toMatrix4()
{
	matrix<T, 4> result;


    T wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;


    // calculate coefficients
    x2 = x + x; y2 = y + y;
    z2 = z + z;
    xx = x * x2; xy = x * y2; xz = x * z2;
    yy = y * y2; yz = y * z2; zz = z * z2;
    wx = w * x2; wy = w * y2; wz = w * z2;

	result(0, 0) = (T)1.0 - (yy + zz);
	result(0, 1) = xy - wz;

	result(0, 2) = xz + wy;

	result(1, 0) = xy + wz;
	result(1, 1) = (T)1.0 - (xx + zz);

	result(1, 2) = yz - wx;

	result(2, 0) = xz - wy;
	result(2, 1) = yz + wx;

	result(2, 2) = (T)1.0 - (xx + yy);

	result(3, 3) = (T)1.0;

	return result;
}

template<typename T>
quaternion<T>::quaternion()
{
	w=0.0; x=0.0; y=0.0; z=0.0;
}

template<typename T>
quaternion<T>::quaternion(T c, T i, T j, T k)
{
	w=c;  x=i;  y=j;  z=k;
}

}	//end rts namespace

//#if __GNUC__ > 3 && __GNUC_MINOR__ > 7
//template<class T> using rtsQuaternion = rts::quaternion<T>;
//#endif

#endif

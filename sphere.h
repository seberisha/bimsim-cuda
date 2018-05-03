#ifndef SPHERESTRUCT_H
#define SPHERESTRUCT_H

#include <ostream>
#include <sstream>
#include <vector>
#include <complex>
#include "fieldslice.h"
#include "dataTypes.h"
#include <algorithm>



struct sphere
{

    //sphere position
    bsPoint p;

    //sphere radius
    ptype a;

    //sphere material index
    unsigned int iMaterial;

    //GPU pointer to the scattered field produced by a plane wave
    //  this is a function of cos(theta) and |r| (distance from sphere center)
    fieldslice Usp;
    fieldslice Uip;
    ptype d_min;
    ptype d_max;

	//sphere order
	int Nl;

	//refractive index for the current lambda
	bsComplex n;

	//external scattering coefficients
	std::vector<bsComplex> B;

	//internal scattering coefficients
	std::vector<bsComplex> A;

    sphere(ptype x = 0.0f, ptype y = 0.0f, ptype z = 0.0f, ptype a = 0.0f, int m = 0, int ang = 128)
    {
        this->p = bsPoint(x, y, z);
        this->a = a;
        this->iMaterial = m;

		//surface = fieldslice(ang, ang/2);
    }

	//assignment operator
	sphere & operator=(const sphere &rhs);

	//copy constructor
	sphere(const sphere &rhs);

	std::string toStr()
	{
		std::stringstream ss;

		ss<<p<<", "<<a<<", "<<iMaterial;

		return ss.str();
	}

	//compute the order required to represent the scattered field

	void calcNl(ptype lambda)
	{
        Nl = ceil( (2 * PI * a) / lambda + 4 * pow( (2 * PI * a) / lambda, (ptype)(1.0/3.0)) + 2);
	}

    //compute the scattering coefficients
	void calcCoeff(ptype lambda, bsComplex n);

	//compute the bessel function look-up tables
	void calcLut(bsComplex* j, bsComplex* h, ptype lambda, bsComplex n, int aR, int rR);
    void calcBesselLut(bsComplex* j, ptype k, bsComplex n, int aR);
	void calcHankelLut(bsComplex* h, ptype k, int rR);

	//calculate the scattering domain Us(theta, r)
	void calcUp(ptype lambda, bsComplex n, rts::quad<ptype, 3> nfPlane, unsigned int R);

	void scalarUsp(bsComplex* h, int rR, int thetaR);
	void scalarUip(bsComplex* j, int aR, int thetaR);



};

#endif

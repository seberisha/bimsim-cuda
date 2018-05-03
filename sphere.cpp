#include "sphere.h"
#include "defaults.h"

#include "rts/math/complex.h"
#include <complex>
#include <stdlib.h>
#include <fstream>

//using namespace rts;
using namespace std;

int cbessjyva(double v,complex<double> z,double &vm,complex<double>*cjv,
    complex<double>*cyv,complex<double>*cjvp,complex<double>*cyvp);

int cbessjyva_sph(int v,complex<double> z,double &vm,complex<double>*cjv,
    complex<double>*cyv,complex<double>*cjvp,complex<double>*cyvp);

int bessjyv_sph(int v, double z, double &vm, double* cjv,
    double* cyv, double* cjvp, double* cyvp);

void sphere::calcCoeff(ptype lambda, bsComplex ri)
{
    /*  These calculations are done at high-precision on the CPU
        since they are only required once for each wavelength.

        Input:

        lambda  =   wavelength of the incident field
        n       =   complex refractive index of the sphere
    */

    //clear the previous coefficients
    A.clear();
    B.clear();

	//convert to an std complex value
	complex<double> nc(ri.real(), ri.imag());
	n = ri;

    //compute the magnitude of the k vector
    double k = 2 * PI / lambda;
    complex<double> kna = nc * k * (double)a;

    //compute the arguments k*a and k*n*a
    complex<double> ka(k * a, 0.0);

    //allocate space for the Bessel functions of the first and second kind (and derivatives)
    int bytes = sizeof(complex<double>) * (Nl + 1);
    complex<double>* cjv_ka = (complex<double>*)malloc(bytes);
    complex<double>* cyv_ka = (complex<double>*)malloc(bytes);
    complex<double>* cjvp_ka = (complex<double>*)malloc(bytes);
    complex<double>* cyvp_ka = (complex<double>*)malloc(bytes);
    complex<double>* cjv_kna = (complex<double>*)malloc(bytes);
    complex<double>* cyv_kna = (complex<double>*)malloc(bytes);
    complex<double>* cjvp_kna = (complex<double>*)malloc(bytes);
    complex<double>* cyvp_kna = (complex<double>*)malloc(bytes);

    //allocate space for the spherical Hankel functions and derivative
    complex<double>* chv_ka = (complex<double>*)malloc(bytes);
    complex<double>* chvp_ka = (complex<double>*)malloc(bytes);

    //compute the bessel functions using the CPU-based algorithm
    double vm;
    cbessjyva_sph(Nl, ka, vm, cjv_ka, cyv_ka, cjvp_ka, cyvp_ka);
    cbessjyva_sph(Nl, kna, vm, cjv_kna, cyv_kna, cjvp_kna, cyvp_kna);

    //compute A for each order
    complex<double> i(0, 1);
    complex<double> a, b, c, d;
    complex<double> An, Bn;
    for(int l=0; l<=Nl; l++)
    {
        //compute the Hankel functions from j and y
        chv_ka[l] = cjv_ka[l] + i * cyv_ka[l];
        chvp_ka[l] = cjvp_ka[l] + i * cyvp_ka[l];

        //Compute A (internal scattering coefficient)
        //compute the numerator and denominator for A
        a = cjv_ka[l] * chvp_ka[l] - cjvp_ka[l] * chv_ka[l];
        b = cjv_kna[l] * chvp_ka[l] - chv_ka[l] * cjvp_kna[l] * nc;

        //calculate A and add it to the list
        An = (2.0 * l + 1.0) * pow(i, l) * (a / b);
        A.push_back(bsComplex(An.real(), An.imag()));


        //Compute B (external scattering coefficient)
        c = cjv_ka[l] * cjvp_kna[l] * nc - cjv_kna[l] * cjvp_ka[l];
        d = cjv_kna[l] * chvp_ka[l] - chv_ka[l] * cjvp_kna[l] * nc;

        //calculate B and add it to the list
        Bn = (2.0 * l + 1.0) * pow(i, l) * (c / d);
        B.push_back(bsComplex(Bn.real(), Bn.imag()));


    }
}

void sphere::calcBesselLut(bsComplex* j, ptype k, bsComplex n, int aR)
{
    /*Compute the look-up-table for spherical bessel functions used inside of the sphere
        j    =   (Nl + 1) x aR array of values
        aR      =   resolution of j
    */

    //allocate space for the Bessel functions of the first and second kind (and derivatives -- which will be ignored)
    int bytes = sizeof(complex<double>) * (Nl + 1);
    complex<double>* cjv_knr = (complex<double>*)malloc(bytes);
    complex<double>* cyv_knr = (complex<double>*)malloc(bytes);
    complex<double>* cjvp_knr = (complex<double>*)malloc(bytes);
    complex<double>* cyvp_knr = (complex<double>*)malloc(bytes);

    //compute the bessel functions using the CPU-based algorithm
    double vm;

    //for each sample along r
    ptype dr = a / (aR - 1);
    ptype r;
    for(int ir = 0; ir < aR; ir++)
    {
        r = ir * dr;
        complex<double> knr( (k*n*r).real(), (k*n*r).imag() );
        cbessjyva_sph(Nl, knr, vm, cjv_knr, cyv_knr, cjvp_knr, cyvp_knr);

        //copy the double data to the bsComplex array
        for(int l=0; l<=Nl; l++)
		{
			//deal with the NaN case at the origin
			if(ir == 0)
			{
				if(l == 0)
					j[ir * (Nl+1)] = 1;
				else
					j[ir * (Nl+1) + l] = 0;
			}
			else
				j[ir * (Nl+1) + l] = bsComplex(cjv_knr[l].real(), cjv_knr[l].imag());
		}
    }

	/*ofstream outfile("besselout.txt");
    for(int ir = 0; ir < aR; ir++)
    {
        for(int l = 0; l<Nl+1; l++)
        {
            outfile<<j[ir * (Nl+1) + l].real()<<"     ";
        }
        outfile<<endl;
    }
	outfile.close();*/

}

void sphere::calcHankelLut(bsComplex* h, ptype k, int rR)
{
	/*Compute the look-up-table for spherical bessel functions used inside of the sphere
        h_out   =   (Nl + 1) x aR array of values
		rmin	=	minimum value of r
		d_max	=	maximum value of r
        rR      =   resolution of h_out
    */

    //allocate space for the Bessel functions of the first and second kind (and derivatives -- which will be ignored)
    int bytes = sizeof(double) * (Nl + 1);
    double* cjv_kr = (double*)malloc(bytes);
    double* cyv_kr = (double*)malloc(bytes);
    double* cjvp_kr = (double*)malloc(bytes);
    double* cyvp_kr = (double*)malloc(bytes);

    //compute the bessel functions using the CPU-based algorithm
    double vm;



    //for each sample along r
    ptype dr = (d_max - max(a, d_min)) / (rR - 1);
    ptype r;
    for(int ir = 0; ir < rR; ir++)
    {
        r = ir * dr + max(a, d_min);
        double kr = k*r;
        bessjyv_sph(Nl, kr, vm, cjv_kr, cyv_kr, cjvp_kr, cyvp_kr);

        //copy the double data to the bsComplex array
        for(int l=0; l<=Nl; l++)
		{
			//h[ir * (Nl+1) + l] = bsComplex(cjv_kr[l].real(), cyv_kr[l].real());
			h[ir * (Nl+1) + l] = bsComplex(cjv_kr[l], cyv_kr[l]);
		}
    }

	/*ofstream outfile("hankelout.txt");
    for(int ir = 0; ir < rR; ir++)
    {
		outfile<<ir*dr + max(a, d_min)<<"     ";
        for(int l = 0; l<=0; l++)
        {
            outfile<<h[ir * (Nl+1) + l].real()<<"     "<<h[ir * (Nl+1) + l].imag()<<"     ";
        }
        outfile<<endl;
    }
	outfile.close();*/
}

void sphere::calcLut(bsComplex* j, bsComplex* h, ptype lambda, bsComplex n, int aR, int rR)
{
    /*Compute the look-up-tables for spherical bessel functions used both inside and outside of the sphere.
        j       =   (Nl + 1) x aR array of values
        j       =   (Nl + 1) x rR array of values
        d_max    =   maximum distance for the LUT
        aR      =   resolution of j_in
        rR      =   resolution of j_out
    */

    //compute the magnitude of the k vector
    double k = 2 * PI / lambda;

	calcBesselLut(j, k, n, aR);
	calcHankelLut(h, k, rR);
}

void sphere::calcUp(ptype lambda, bsComplex n, rts::quad<ptype, 3> nfPlane, unsigned int R)
{
    //calculate the parameters of the lookup table

    //first find the distance to the closest and furthest points on the nearfield plane
    d_min = nfPlane.dist(p);
    d_max = nfPlane.dist_max(p);

    //compute the radius of the cross-section of the sphere with the plane
    ptype a_inter = 0;
    if(d_min < a)
        a_inter = sqrt(a - d_min);


	//calculate the resolution of the Usp and Uip lookup tables
	int aR = 1 + 2 * R * a_inter / (nfPlane(0, 0) - nfPlane(1, 1)).len();
	int dR = 2 * R;
	int thetaR = DEFAULT_SPHERE_THETA_R;

	//allocate space for the bessel function LUTs
	bsComplex* j = (bsComplex*)malloc(sizeof(bsComplex) * (Nl + 1) * aR);
	bsComplex* h = (bsComplex*)malloc(sizeof(bsComplex) * (Nl + 1) * dR);

	calcLut(j, h, lambda, n, aR, dR);

	//allocate space for the Usp lookup texture
	Usp.R[0] = dR;
	Usp.R[1] = thetaR;
	Usp.init_gpu();

	//allocate space for the Uip lookup texture
	Uip.R[0] = aR;
	Uip.R[1] = thetaR;
	Uip.init_gpu();



	scalarUsp(h, dR, thetaR);
	scalarUip(j, aR, thetaR);

	scalarslice UspMag = Usp.Mag();
    //UspMag.toImage("Usp.bmp", true);

	scalarslice UipMag = Uip.Mag();
    //UipMag.toImage("Uip.bmp", true);

	//free memory
	free(j);
	free(h);

}

sphere& sphere::operator=(const sphere &rhs)
{
	p = rhs.p;
	a = rhs.a;
	iMaterial = rhs.iMaterial;
	Nl = rhs.Nl;
	n = rhs.n;
	B = rhs.B;
	A = rhs.A;

	return *this;
}

sphere::sphere(const sphere &rhs)
{
	p = rhs.p;
	a = rhs.a;
	iMaterial = rhs.iMaterial;
	Nl = rhs.Nl;
	n = rhs.n;
	B = rhs.B;
	A = rhs.A;
}

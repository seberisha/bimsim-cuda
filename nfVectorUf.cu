#include "nearfield.h"
#include "rts/math/spherical_bessel.h"
#include "rts/math/legendre.h"
#include <stdlib.h>
#include "rts/cuda/error.h"
#include "rts/cuda/timer.h"

//Incident field for a single plane wave
__global__ void gpuVectorUfp(bsComplex* Uf, bsVector k, ptype kmag, bsPoint f, ptype A, bsRect ABCD, int uR, int vR)
{
	/*Compute the scalar focused field using Debye focusing
		k		= direction of focused light, where |k| = 2*pi/lambda
		P		= rect struct describing the field slice
		rX, rY	= resolution of the field slice
		cNAin	= inner NA of the condenser
		cNAout	= outer NA of the condenser
	*/

	//get the current coordinate in the plane slice
	int iu = blockIdx.x * blockDim.x + threadIdx.x;
	int iv = blockIdx.y * blockDim.y + threadIdx.y;

	//make sure that the thread indices are in-bounds
	if(iu >= uR || iv >= vR) return;

	//compute the index (easier access to the scalar field array)
	int i = iv*uR + iu;

	//compute the parameters for u and v
	ptype u = (ptype)iu / uR;
	ptype v = (ptype)iv / vR;

	//get the rtsPoint in world space and then the r vector
	bsPoint p = ABCD(u, v);
	bsVector r = p - f;
	//ptype d = r.len();

	ptype k_dot_r = kmag * k.dot(r);
	bsComplex d(0, k_dot_r);

	Uf[i] = exp(d) * A;

}

//Incident field for a focused point source
__global__ void gpuVectorUf(bsComplex* Uf, bsVector k, ptype kmag, bsPoint f, ptype A, bsRect ABCD, int uR, int vR, ptype cosAlpha, ptype cosBeta, int nl, ptype j_conv = 1.4)
{
    //Compute the scalar focused field using Debye focusing
	//	k		= direction of focused light, where |k| = 2*pi/lambda
	//	P		= rect struct describing the field slice
	//	rX, rY	= resolution of the field slice
	//	cNAin	= inner NA of the condenser
	//	cNAout	= outer NA of the condenser


	//get the current coordinate in the plane slice
	int iu = blockIdx.x * blockDim.x + threadIdx.x;
	int iv = blockIdx.y * blockDim.y + threadIdx.y;

	//make sure that the thread indices are in-bounds
	if(iu >= uR || iv >= vR) return;

	//compute the index (easier access to the scalar field array)
	int i = iv*uR + iu;

	//compute the parameters for u and v
	ptype u = (ptype)iu / (uR);
	ptype v = (ptype)iv / (vR);

	//get the rtsPoint in world space and then the r vector
	bsPoint p = ABCD(u, v);
	bsVector r = p - f;
	ptype d = r.len();
	if(d < EPSILON_FLOAT)
	{
        Uf[i] = A * 2 * PI * (cosAlpha - cosBeta);
        return;
    }

	//get info for the light direction and frequency
	//k = k.norm();
	r = r.norm();

	//compute the imaginary factor i^l
	bsComplex im = bsComplex(0, 1);
	bsComplex il = bsComplex(1, 0);

	//Bessel and Legendre functions are computed dynamically to save memory
	//initialize the Bessel and Legendre functions
	ptype j[2];
	ptype kd = kmag * d;
	rts::init_sbesselj<ptype>(kd, j);

	ptype P[2];
	//get the angle between k and r (light direction and position vector)
	ptype cosTheta;
	cosTheta = k.dot(r);

	//deal with the degenerate case where r == 0
	//if(isnan(cosTheta))
    //    cosTheta = 0;
	rts::init_legendre<ptype>(cosTheta, P[0], P[1]);

	//initialize legendre functions for the cassegrain angles
	ptype Palpha[3];
	//ptype cosAlpha = cos(asin(cNAin));
	rts::init_legendre<ptype>(cosAlpha, Palpha[0], Palpha[1]);
	Palpha[2] = 1;

	ptype Pbeta[3];
	//ptype cosBeta = cos(asin(cNAout));
	rts::init_legendre<ptype>(cosBeta, Pbeta[0], Pbeta[1]);
	Pbeta[2] = 1;

	//for each order l
	bsComplex sumUf(0.0, 0.0);
	ptype jl = 0.0;
	ptype Pl;
	for(int l = 0; l<=nl; l++)
	{

		if(l==0)
		{

			jl = j[0];
			Pl = P[0];
		}
		else if(l==1)
		{
			jl = j[1];
			Pl = P[1];

			//adjust the cassegrain Legendre function
			Palpha[2] = Palpha[0];
			rts::shift_legendre<ptype>(l+1, cosAlpha, Palpha[0], Palpha[1]);
			Pbeta[2] = Pbeta[0];
			rts::shift_legendre<ptype>(l+1, cosBeta, Pbeta[0], Pbeta[1]);
		}
		else
		{
			rts::shift_sbesselj<ptype>(l, kd, j);//, j_conv);
			rts::shift_legendre<ptype>(l, cosTheta, P[0], P[1]);

			jl = j[1];
			Pl = P[1];

			//adjust the cassegrain outer Legendre function
			Palpha[2] = Palpha[0];
			rts::shift_legendre<ptype>(l+1, cosAlpha, Palpha[0], Palpha[1]);
			Pbeta[2] = Pbeta[0];
			rts::shift_legendre<ptype>(l+1, cosBeta, Pbeta[0], Pbeta[1]);
		}

		sumUf += il * jl * Pl * (Palpha[1] - Palpha[2] - Pbeta[1] + Pbeta[2]);

		il *= im;
	}

	Uf[i] = sumUf * 2 * PI * A;

}


void nearfieldStruct::vectorUf()
{


    gpuStartTimer();

	//create one thread for each pixel of the field slice
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
	dim3 dimGrid((Uf.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (Uf.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

	//if we are computing a plane wave, call the gpuScalarUfp function
	if(planeWave)
	{
      //  std::cout<<"Calculating vector plane wave..."<<std::endl;
		gpuVectorUfp<<<dimGrid, dimBlock>>>(Uf.x_hat, k, 2 * PI / lambda, focus, A, pos, Uf.R[0], Uf.R[1]);
	}
	//otherwise compute the condenser info and create a focused field
	else
	{
		//pre-compute the cosine of the obscuration and objective angles
		ptype cosAlpha = cos(asin(condenser[0]));
		ptype cosBeta = cos(asin(condenser[1]));
		//compute the scalar Uf field (this will be in the x_hat channel of Uf)
		gpuVectorUf<<<dimGrid, dimBlock>>>(Uf.x_hat, k, 2 * PI / lambda, focus, A, pos, Uf.R[0], Uf.R[1], cosAlpha, cosBeta, m);
	}

	t_Uf = gpuStopTimer();
}

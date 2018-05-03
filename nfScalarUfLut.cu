#include "nearfield.h"

#include "rts/math/legendre.h"
#include "rts/cuda/error.h"
#include "rts/cuda/timer.h"

texture<float, cudaTextureType2D> texJ;

__global__ void gpuScalarUfp(bsComplex* Uf, bsVector k, ptype kmag, bsPoint f, ptype A, bsRect ABCD, int uR, int vR);

__global__ void gpuScalarUfLut(bsComplex* Uf, bsRect ABCD, int uR, int vR, bsPoint f, bsVector k, ptype A, ptype cosAlpha, ptype cosBeta, int nl, ptype dmin, ptype dmax, int dR)
{
    /*This function computes the focused field for a 2D slice

    Uf      =   destination field slice
    ABCD    =   plane representing the field slice in world space
    uR, vR  =   resolution of the Uf field
    f       =   focal point of the condenser
    k       =   direction of the incident light
    A       =   amplitude of the incident field
    cosAlpha=   cosine of the solid angle subtended by the condenser obscuration
    cosBeta =   cosine of the solid angle subtended by the condenser aperature
    nl      =   number of orders used to compute the field
    dR      =   number of Bessel function values in the look-up texture

    */

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

	if(d == 0)
	{
        Uf[i] = A * 2 * PI * (cosAlpha - cosBeta);
        return;
    }

	//get info for the light direction and frequency
	r = r.norm();

	//compute the imaginary factor i^l
	bsComplex im = bsComplex(0, 1);
	bsComplex il = bsComplex(1, 0);

	//Legendre functions are computed dynamically to save memory
	//initialize the Legendre functions

	ptype P[2];
	//get the angle between k and r (light direction and position vector)
	ptype cosTheta;
	cosTheta = k.dot(r);

	rts::init_legendre<ptype>(cosTheta, P[0], P[1]);

	//initialize legendre functions for the cassegrain angles
	ptype Palpha[3];
	rts::init_legendre<ptype>(cosAlpha, Palpha[0], Palpha[1]);
	Palpha[2] = 1;

	ptype Pbeta[3];
	rts::init_legendre<ptype>(cosBeta, Pbeta[0], Pbeta[1]);
	Pbeta[2] = 1;

	//for each order l
	bsComplex sumUf(0, 0);
	ptype jl = 0;
	ptype Pl;
	ptype di = ( (d - dmin)/(dmax - dmin) ) * (dR - 1);
	for(int l = 0; l<=nl; l++)
	{
        jl = tex2D(texJ, l + 0.5f, di + 0.5f);
		if(l==0)
			Pl = P[0];
		else if(l==1)
		{
			Pl = P[1];

			//adjust the cassegrain Legendre function
			Palpha[2] = Palpha[0];
			rts::shift_legendre<ptype>(l+1, cosAlpha, Palpha[0], Palpha[1]);
			Pbeta[2] = Pbeta[0];
			rts::shift_legendre<ptype>(l+1, cosBeta, Pbeta[0], Pbeta[1]);
		}
		else
		{
			rts::shift_legendre<ptype>(l, cosTheta, P[0], P[1]);

			Pl = P[1];

			//adjust the cassegrain outer Legendre function
			Palpha[2] = Palpha[0];
			rts::shift_legendre<ptype>(l+1, cosAlpha, Palpha[0], Palpha[1]);
			Pbeta[2] = Pbeta[0];
			rts::shift_legendre<ptype>(l+1, cosBeta, Pbeta[0], Pbeta[1]);
		}

		sumUf += il * jl * Pl * (Palpha[1] - Palpha[2] - Pbeta[1] + Pbeta[2]);
		//sumUf += jl;

		il *= im;
	}

	Uf[i] = sumUf * 2 * PI * A;
	//Uf[i] = u;
	//return;
}

void nearfieldStruct::scalarUfLut()
{
    gpuStartTimer();
	
    //calculate the minimum and maximum points in the focused field
    d_min = pos.dist(focus);
    d_max = pos.dist_max(focus);

    //allocate space for the Bessel function
    int dR = 2 * max(Uf.R[0], Uf.R[1]);
    ptype* j = NULL;
	j = (ptype*) malloc(sizeof(ptype) * dR * (m+1));

	//calculate Bessel function LUT
	calcBesselLut(j, d_min, d_max, dR);
	
    //create a CUDA array structure and specify the format description
	cudaArray* arrayJ;
    cudaChannelFormatDesc channelDesc =
        cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	
    //allocate memory
    HANDLE_ERROR(cudaMallocArray(&arrayJ, &channelDesc, m+1, dR));
	
    //specify texture properties
    texJ.addressMode[0] = cudaAddressModeMirror;
    texJ.addressMode[1] = cudaAddressModeMirror;
    texJ.filterMode     = cudaFilterModeLinear;
    texJ.normalized     = false;

    //bind the texture to the array
    HANDLE_ERROR(cudaBindTextureToArray(texJ, arrayJ, channelDesc));

    //copy the CPU Bessel LUT to the GPU-based array
    HANDLE_ERROR( cudaMemcpy2DToArray(arrayJ, 0, 0, j, (m+1)*sizeof(float), (m+1)*sizeof(float), dR, cudaMemcpyHostToDevice));

    //----------------Compute the focused field
    //create one thread for each pixel of the field slice
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
	dim3 dimGrid((Uf.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (Uf.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

	//if we are computing a plane wave, call the gpuScalarUfp function
	if(planeWave)
	{
		gpuScalarUfp<<<dimGrid, dimBlock>>>(Uf.x_hat, k, 2 * PI / lambda, focus, A, pos, Uf.R[0], Uf.R[1]);
	}
	//otherwise compute the condenser info and create a focused field
	else
	{
		//pre-compute the cosine of the obscuration and objective angles
		ptype cosAlpha = cos(asin(condenser[0]));
		ptype cosBeta = cos(asin(condenser[1]));
		//compute the scalar Uf field (this will be in the x_hat channel of Uf)
		gpuScalarUfLut<<<dimGrid, dimBlock>>>(Uf.x_hat, pos, Uf.R[0], Uf.R[1], focus, k, A, cosAlpha, cosBeta, m, d_min, d_max, dR);
	}

	
    //free everything
	free(j);
	
	HANDLE_ERROR(cudaFreeArray(arrayJ));

	t_Uf = gpuStopTimer();
}

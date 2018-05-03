#include "sphere.h"
#include "rts/math/legendre.h"

__global__ void gpuScalarUsp(bsComplex* Usp, bsComplex* h, bsComplex* B, int Nl, int rR, int thetaR)
{
    //get the current coordinate in the plane slice
	int ir = blockIdx.x * blockDim.x + threadIdx.x;
	int itheta = blockIdx.y * blockDim.y + threadIdx.y;

	//make sure that the thread indices are in-bounds
	if(itheta >= thetaR || ir >= rR) return;

	int i = itheta * rR + ir;

	//ptype dr = (rmax - a) / (rR - 1);
	ptype dtheta = (PI) / (thetaR - 1);

	//comptue the current angle and distance
	//ptype r = dr * ir + a;
	ptype theta = dtheta * itheta;
	ptype cos_theta = cos(theta);

	//initialize the Legendre polynomial
	ptype P[2];
	rts::init_legendre<ptype>(cos_theta, P[0], P[1]);

	//initialize the result
	bsComplex Us((ptype)0, (ptype)0);

    //for each order l
    for(int l=0; l <= Nl; l++)
    {
        if(l == 0)
        {
            Us += B[l] * h[ir * (Nl+1) + l] * P[0];
            //Us += P[0];
        }
        else
        {
            if(l > 1)
            {
                rts::shift_legendre<ptype>(l, cos_theta, P[0], P[1]);
            }
            Us += B[l] * h[ir * (Nl+1) + l] * P[1];
            //Us += P[1];
        }


    }
	Usp[i] = Us;
	//Usp[i] = h[ir * (Nl+1)];
	//Usp[i] = ir;

}

__global__ void gpuScalarUip(bsComplex* Uip, bsComplex* j, bsComplex* A, int Nl, int aR, int thetaR)
{
    //get the current coordinate in the plane slice
	int ia = blockIdx.x * blockDim.x + threadIdx.x;
	int itheta = blockIdx.y * blockDim.y + threadIdx.y;

	//make sure that the thread indices are in-bounds
	if(itheta >= thetaR || ia >= aR) return;

	int i = itheta * aR + ia;

	ptype dtheta = (PI) / (thetaR - 1);

	//comptue the current angle and distance
	ptype theta = dtheta * itheta;
	ptype cos_theta = cos(theta);

	//initialize the Legendre polynomial
	ptype P[2];
	rts::init_legendre<ptype>(cos_theta, P[0], P[1]);

	//initialize the result
	bsComplex Ui((ptype)0, (ptype)0);

    //for each order l
    for(int l=0; l <= Nl; l++)
    {
        if(l == 0)
        {
            Ui += A[l] * j[ia * (Nl+1) + l] * P[0];
        }
        else
        {
            if(l > 1)
            {
                rts::shift_legendre<ptype>(l, cos_theta, P[0], P[1]);
            }
            Ui += A[l] * j[ia * (Nl+1) + l] * P[1];
        }


    }
	Uip[i] = Ui;
}

void sphere::scalarUsp(bsComplex* h, int rR, int thetaR)
{
	//copy the hankel function to the GPU
    bsComplex* gpu_h;
    HANDLE_ERROR( cudaMalloc( (void**)&gpu_h, sizeof(bsComplex) * (Nl + 1) * rR ) );
    HANDLE_ERROR( cudaMemcpy( gpu_h, h, sizeof(bsComplex) * (Nl + 1) * rR, cudaMemcpyHostToDevice ) );

    //allocate memory for the scattering coefficients
    bsComplex* gpuB;
    HANDLE_ERROR(cudaMalloc((void**) &gpuB, (Nl+1) * sizeof(bsComplex)));
    //copy the scattering coefficients to the GPU
    HANDLE_ERROR(cudaMemcpy(gpuB, &B[0], (Nl+1) * sizeof(bsComplex), cudaMemcpyHostToDevice));

    //create one thread for each pixel of the field slice
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
	dim3 dimGrid((Usp.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (Usp.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

	gpuScalarUsp<<<dimGrid, dimBlock>>>(Usp.x_hat, gpu_h, gpuB, Nl, rR, thetaR);

	//free memory
	cudaFree(gpu_h);
	cudaFree(gpuB);

}

void sphere::scalarUip(bsComplex* j, int rR, int thetaR)
{
	//copy the bessel and hankel LUTs to the GPU
    bsComplex* gpu_j;
    HANDLE_ERROR( cudaMalloc( (void**)&gpu_j, sizeof(bsComplex) * (Nl + 1) * rR ) );
    HANDLE_ERROR( cudaMemcpy( gpu_j, j, sizeof(bsComplex) * (Nl + 1) * rR, cudaMemcpyHostToDevice ) );

    //allocate memory for the scattering coefficients
    bsComplex* gpuA;
    HANDLE_ERROR(cudaMalloc((void**) &gpuA, (Nl+1) * sizeof(bsComplex)));
    //copy the scattering coefficients to the GPU
    HANDLE_ERROR(cudaMemcpy(gpuA, &A[0], (Nl+1) * sizeof(bsComplex), cudaMemcpyHostToDevice));

    //create one thread for each pixel of the field slice
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
	dim3 dimGrid((Uip.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (Uip.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

	gpuScalarUip<<<dimGrid, dimBlock>>>(Uip.x_hat, gpu_j, gpuA, Nl, rR, thetaR);

	//free memory
	cudaFree(gpu_j);
	cudaFree(gpuA);

}

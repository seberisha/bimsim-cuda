#include "nearfield.h"
#include <stdlib.h>
#include "rts/cuda/error.h"
#include "rts/cuda/timer.h"

__global__ void gpuScalarUsp(bsComplex* Ufx, bsComplex* Ufy, bsComplex* Ufz,
							 bsComplex* Ux, bsComplex* Uy, bsComplex* Uz,
							 bsPoint* ps, ptype* as, int ns, bsRect ABCD, int uR, int vR)
{

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
	bsVector r;
	ptype d;

	//if we are inside of a sphere, return
	for(int is=0; is<ns; is++)
	{
		r = p - ps[is];
		d = r.len();
		if(d < as[is])
		{
			//printf("\t\t sumUf inside sphere -- d: %f\n",d);
			return;
		}
	}

	//otherwise add the focused field to the full field
	if(Ufx != NULL)
		Ux[i] += Ufx[i];
	if(Ufy != NULL)
		Uy[i] += Ufy[i];
	if(Ufz != NULL)
		Uz[i] += Ufz[i];
}

void nearfieldStruct::sumUf()
{


	//create arrays to store sphere positions and radii
	int nSpheres = sVector.size();

	//if the number of spheres is zero, just copy the incident field
	if(nSpheres == 0)
	{
		if(U.x_hat != NULL)
			HANDLE_ERROR(cudaMemcpy(U.x_hat, Uf.x_hat, sizeof(bsComplex) * U.R[0] * U.R[1], cudaMemcpyDeviceToDevice));
		if(U.y_hat != NULL)
			HANDLE_ERROR(cudaMemcpy(U.y_hat, Uf.y_hat, sizeof(bsComplex) * U.R[0] * U.R[1], cudaMemcpyDeviceToDevice));
		if(U.z_hat != NULL)
			HANDLE_ERROR(cudaMemcpy(U.z_hat, Uf.z_hat, sizeof(bsComplex) * U.R[0] * U.R[1], cudaMemcpyDeviceToDevice));
		return;
	}

	//time the calculation of the focused field
	//gpuStartTimer();

	bsPoint* cpu_p = (bsPoint*)malloc(sizeof(bsPoint) * nSpheres);
	ptype* cpu_a = (ptype*)malloc(sizeof(ptype) * nSpheres);

	//copy the sphere positions and radii to the new arrays
	for(int s=0; s<nSpheres; s++)
	{
		cpu_p[s] = sVector[s].p;
		cpu_a[s] = sVector[s].a;
	}

	//copy the arrays to the gpu
	bsPoint* gpu_p;
	HANDLE_ERROR(cudaMalloc( (void**) &gpu_p, sizeof(bsPoint) * nSpheres));
	HANDLE_ERROR(cudaMemcpy(gpu_p, cpu_p, sizeof(bsPoint) * nSpheres, cudaMemcpyHostToDevice));
	ptype* gpu_a;
	HANDLE_ERROR(cudaMalloc( (void**) &gpu_a, sizeof(ptype) * nSpheres));
	HANDLE_ERROR(cudaMemcpy(gpu_a, cpu_a, sizeof(ptype) * nSpheres, cudaMemcpyHostToDevice));


	//create one thread for each pixel of the field slice
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
	dim3 dimGrid((U.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (U.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

	//copy the focused field
	gpuScalarUsp<<<dimGrid, dimBlock>>>(Uf.x_hat,
										Uf.y_hat,
										Uf.z_hat,
										U.x_hat,
										U.y_hat,
										U.z_hat,
										gpu_p,
										gpu_a,
										nSpheres,
										pos,
										U.R[0],
										U.R[1]);



    //free sphere lists
    HANDLE_ERROR(cudaFree(gpu_p));
    HANDLE_ERROR(cudaFree(gpu_a));


}

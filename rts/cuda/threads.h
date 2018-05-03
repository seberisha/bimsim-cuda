#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_callable.h"

#ifndef CUDA_THREADS_H
#define CUDA_THREADS_H

#define MAX_GRID        65535

__device__ unsigned int ThreadIndex1D()
{
    return blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
}

dim3 GenGrid1D(unsigned int N, unsigned int blocksize = 128)
{
    dim3 dimgrid;

    dimgrid.x = (N + blocksize - 1)/blocksize;
    dimgrid.y = 1;
    dimgrid.z = 1;

    if(dimgrid.x > MAX_GRID)
    {
        dimgrid.y = (dimgrid.x + MAX_GRID - 1) / MAX_GRID;
        dimgrid.x = MAX_GRID;
    }

    return dimgrid;

}


#endif

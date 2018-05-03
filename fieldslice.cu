#include "fieldslice.h"
#include "dataTypes.h"
#include "rts/cuda/error.h"
#include "rts/cuda/threads.h"

__global__ void field_intensity(bsComplex* x, bsComplex* y, bsComplex* z, ptype* I, unsigned int N)
{
    //compute the index for this thread
	//int i = blockIdx.x * blockDim.x + threadIdx.x;
	int i = ThreadIndex1D();

	if(i >= N) return;

	ptype xm = x[i].abs();


	if(y != NULL && z != NULL)
	{
		ptype ym = y[i].abs();
		ptype zm = z[i].abs();
		I[i] = xm*xm + ym*ym + zm*zm;
	}
	else
	{
		I[i] = xm*xm;
	}
}


__global__ void resample_intensity_without_integrating(bsComplex* x, bsComplex* y, bsComplex* z, ptype* D, int uR, int vR, int ss)
{
    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index into the detector
    int i = iv*uR + iu;

    //compute the index into the field
    int fi;

    //initialize the intensity for the pixel to zero
    ptype I = 0;
    ptype xm = 0;
    ptype ym = 0;
    ptype zm = 0;

    int ix, iy;
    for(ix = 0; ix<ss; ix++)
        for(iy = 0; iy<ss; iy++)
        {
            //fi = iv*ss*ss*uR + iy*ss*uR + iu*ss + ix;
            fi = (iv*ss + iy)*ss*uR + iu*ss + ix;
            if(x !=NULL)
                xm = x[fi].abs();
            if(y != NULL)
                ym = y[fi].abs();
            if(z != NULL)
                zm = z[fi].abs();
            I += xm*xm + ym*ym + zm*zm;
        }

    D[i] = I/(ss*ss);
   // if (i<128)
}


__global__ void resample_intensity(bsComplex* x, bsComplex* y, bsComplex* z, ptype* D, int uR, int vR, int ss, float scale)
{
	//get the current coordinate in the plane slice
	int iu = blockIdx.x * blockDim.x + threadIdx.x;
	int iv = blockIdx.y * blockDim.y + threadIdx.y;

	//make sure that the thread indices are in-bounds
	if(iu >= uR || iv >= vR) return;

	//compute the index into the detector
	int i = iv*uR + iu;

	//compute the index into the field
	int fi;

	//initialize the intensity for the pixel to zero
	ptype I = 0;
	ptype xm = 0;
	ptype ym = 0;
	ptype zm = 0;

	int ix, iy;
	for(ix = 0; ix<ss; ix++)
		for(iy = 0; iy<ss; iy++)
		{
			//fi = iv*ss*ss*uR + iy*ss*uR + iu*ss + ix;
			fi = (iv*ss + iy)*ss*uR + iu*ss + ix;
			if(x !=NULL)
				xm = x[fi].abs();
			if(y != NULL)
				ym = y[fi].abs();
			if(z != NULL)
				zm = z[fi].abs();
			I += xm*xm + ym*ym + zm*zm;
		}

    D[i] += scale*I/(ss*ss);
}

__global__ void field_real(bsComplex* field_component, ptype* V, unsigned int N)
{
    //compute the index for this thread
	int i = ThreadIndex1D();
	if(i >= N) return;

	V[i] = field_component[i].real();
}

__global__ void field_imaginary(bsComplex* field_component, ptype* V, unsigned int N)
{
    //compute the index for this thread
	//int i = blockIdx.x * blockDim.x + threadIdx.x;
	int i = ThreadIndex1D();
	if(i >= N) return;

	V[i] = field_component[i].imag();
}

__global__ void field_sqrt(ptype* input, ptype* output, unsigned int N)
{
	//compute the index for this thread
	//int i = blockIdx.x * blockDim.x + threadIdx.x;
	int i = ThreadIndex1D();
	if(i >= N) return;

	output[i] = sqrt(input[i]);

}


__global__ void field_scale(bsComplex* x, bsComplex* y, bsComplex* z, unsigned int N, ptype v)
{
    //compute the index for this thread
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i >= N) return;

	if(x != NULL)
        x[i] *= v;
    if(y != NULL)
        y[i] *= v;
    if(z != NULL)
        z[i] *= v;
}


scalarslice fieldslice::Mag()
{
	//compute the magnitude of the field at each rtsPoint in the slice

    scalarslice* result = new scalarslice(R[0], R[1]);

	//compute the total number of values in the slice
	unsigned int N = R[0] * R[1];
	//int gridDim = (N+BLOCK-1)/BLOCK;
	dim3 gridDim = GenGrid1D(N, BLOCK);

	field_intensity<<<gridDim, BLOCK>>>(x_hat, y_hat, z_hat, result->S, N);
	field_sqrt<<<gridDim, BLOCK>>>(result->S, result->S, N);

	return *result;
}

scalarslice fieldslice::Real()
{
	//compute the magnitude of the field at each rtsPoint in the slice

	//create a scalar slice at the same resolution as the field
	scalarslice* result = new scalarslice(R[0], R[1]);

	//compute the total number of values in the slice
	unsigned int N = R[0] * R[1];
	//int gridDim = (N+BLOCK-1)/BLOCK;
	dim3 gridDim = GenGrid1D(N, BLOCK);

	field_real<<<gridDim, BLOCK>>>(x_hat, result->S, N);

	return *result;
}

scalarslice fieldslice::Imag()
{
	//compute the magnitude of the field at each rtsPoint in the slice

	//create a scalar slice at the same resolution as the field
	scalarslice* result = new scalarslice(R[0], R[1]);

	//compute the total number of values in the slice
	unsigned int N = R[0] * R[1];
	//int gridDim = (N+BLOCK-1)/BLOCK;
	dim3 gridDim = GenGrid1D(N, BLOCK);

	field_imaginary<<<gridDim, BLOCK>>>(x_hat, result->S, N);

	return *result;
}

void fieldslice::IntegrateAndResample(scalarslice* detector, unsigned int supersample, float scale)
{
    //compute the intensity and resample at the detector resolution
	unsigned int D[2];
	D[0] = detector->R[0];
	D[1] = detector->R[1];

	//create one thread for each detector pixel
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
	dim3 dimGrid((D[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (D[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

    resample_intensity<<<dimGrid, dimBlock>>>(x_hat, y_hat, z_hat, detector->S, D[0], D[1], supersample, scale);
}

void fieldslice::ResampleAndSave(scalarslice* detector, unsigned int supersample)
{
    //compute the intensity and resample at the detector resolution
    unsigned int D[2];
    D[0] = detector->R[0];
    D[1] = detector->R[1];

    //create one thread for each detector pixel
    dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    dim3 dimGrid((D[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (D[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

    resample_intensity_without_integrating<<<dimGrid, dimBlock>>>(x_hat, y_hat, z_hat, detector->S, D[0], D[1], supersample);
}


scalarslice fieldslice::Intensity()
{
	//compute the magnitude of the field at each rtsPoint in the slice

	//create a scalar slice at the same resolution as the field
	scalarslice* result = new scalarslice(R[0], R[1]);

	//compute the total number of values in the slice
	unsigned int N = R[0] * R[1];
	int gridDim = (N+BLOCK-1)/BLOCK;

	field_intensity<<<gridDim, BLOCK>>>(x_hat, y_hat, z_hat, result->S, N);

	return *result;
}

void fieldslice::ScaleField(ptype v)
{
    //This function scales the field by some constant value v
    //This is mostly used for the inverse FFT, which has to divide the field by R^2

    //compute the total number of values in the slice
	unsigned int N = R[0] * R[1];
	int gridDim = (N+BLOCK-1)/BLOCK;

	field_scale<<<gridDim, BLOCK>>>(x_hat, y_hat, z_hat, N, v);

}

void fieldslice::init_gpu()
{
	//if the field has no size, return
	if(R[0] == 0 || R[1] == 0)
		return;

    //free any previous memory allocations
    if(x_hat)
        HANDLE_ERROR(cudaFree(x_hat));
    if(y_hat)
        HANDLE_ERROR(cudaFree(y_hat));
    if(z_hat)
        HANDLE_ERROR(cudaFree(z_hat));

    //allocate space on the GPU for the field slice
	HANDLE_ERROR(cudaMalloc((void**)&x_hat, R[0] * R[1] * sizeof(bsComplex)));

	if(!scalarField)
	{
		HANDLE_ERROR(cudaMalloc((void**)&y_hat, R[0] * R[1] * sizeof(bsComplex)));
		//HANDLE_ERROR(cudaMemset(y_hat, 0, R[0] * R[1] * sizeof(bsComplex)));

		HANDLE_ERROR(cudaMalloc((void**)&z_hat, R[0] * R[1] * sizeof(bsComplex)));
		//HANDLE_ERROR(cudaMemset(z_hat, 0, R[0] * R[1] * sizeof(bsComplex)));
	}

	clear_gpu();
}

void fieldslice::kill_gpu()
{
    if(x_hat != NULL)
        HANDLE_ERROR(cudaFree(x_hat));
    if(y_hat != NULL)
        HANDLE_ERROR(cudaFree(y_hat));
    if(z_hat != NULL)
        HANDLE_ERROR(cudaFree(z_hat));

	x_hat = y_hat = z_hat = NULL;

}

void fieldslice::clear_gpu()
{
	int memsize = R[0] * R[1] * sizeof(bsComplex);
	if(x_hat != NULL)
		HANDLE_ERROR(cudaMemset(x_hat, 0, memsize));
	if(y_hat != NULL)
		HANDLE_ERROR(cudaMemset(y_hat, 0, memsize));
	if(z_hat != NULL)
		HANDLE_ERROR(cudaMemset(z_hat, 0, memsize));

}

__global__ void copy_crop(bsComplex* source, bsComplex* dest, int u, int v, int su, int sv, int uR, int vR)
{
    //get the current coordinate in the plane slice
	int iu = blockIdx.x * blockDim.x + threadIdx.x;
	int iv = blockIdx.y * blockDim.y + threadIdx.y;

	//make sure that the thread indices are in-bounds
	if(iu >= su || iv >= sv) return;

	//compute the destination index
	int i = iv*su + iu;

	//compute the source index
	int sourceV = v + iv;
	int sourceU = u + iu;
	int is = sourceV * uR + sourceU;

	dest[i] = source[is];

}

fieldslice fieldslice::crop(int u, int v, int su, int sv)
{
	//create a new field slice with the appropriate settings
	fieldslice result(su, sv);
	result.scalarField = scalarField;

	//allocate space for the new field
	//result.init_gpu();

	//create one thread for each pixel of the field slice
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
	dim3 dimGrid((su + SQRT_BLOCK -1)/SQRT_BLOCK, (sv + SQRT_BLOCK - 1)/SQRT_BLOCK);

	//call a kernel to copy the cropped to the new field slice
	if(x_hat != NULL)
		copy_crop<<<dimGrid, dimBlock>>>(x_hat, result.x_hat, u, v, su, sv, R[0], R[1]);
	if(y_hat != NULL)
		copy_crop<<<dimGrid, dimBlock>>>(y_hat, result.y_hat, u, v, su, sv, R[0], R[1]);
	if(z_hat != NULL)
		copy_crop<<<dimGrid, dimBlock>>>(z_hat, result.z_hat, u, v, su, sv, R[0], R[1]);

	return result;
}

fieldslice::fieldslice(const fieldslice& rhs)
{
	R[0] = rhs.R[0];
	R[1] = rhs.R[1];
	scalarField = rhs.scalarField;

	x_hat = y_hat = z_hat = NULL;

	unsigned int bytes = sizeof(bsComplex) * R[0] * R[1];
	if(rhs.x_hat != NULL)
	{
		HANDLE_ERROR(cudaMalloc( (void**)&x_hat, bytes));
		HANDLE_ERROR(cudaMemcpy( x_hat, rhs.x_hat, bytes, cudaMemcpyDeviceToDevice));
	}
	if(rhs.y_hat != NULL)
	{
		HANDLE_ERROR(cudaMalloc( (void**)&y_hat, bytes));
		HANDLE_ERROR(cudaMemcpy( y_hat, rhs.y_hat, bytes, cudaMemcpyDeviceToDevice));
	}
	if(rhs.z_hat != NULL)
	{
		HANDLE_ERROR(cudaMalloc( (void**)&z_hat, bytes));
		HANDLE_ERROR(cudaMemcpy( z_hat, rhs.z_hat, bytes, cudaMemcpyDeviceToDevice));
	}

}

fieldslice& fieldslice::operator=(const fieldslice& rhs)
{
	//make sure this isn't a self-allocation
	if(this != &rhs)
	{
		//make a shallow copy
		R[0] = rhs.R[0];
		R[1] = rhs.R[1];
		scalarField = rhs.scalarField;

		//initialize to new parameters
		init_gpu();

		//make a deep copy
		unsigned int bytes = sizeof(bsComplex) * R[0] * R[1];
		if(x_hat != NULL)
			HANDLE_ERROR(cudaMemcpy(x_hat, rhs.x_hat, bytes, cudaMemcpyDeviceToDevice));
		if(y_hat != NULL)
			HANDLE_ERROR(cudaMemcpy(y_hat, rhs.y_hat, bytes, cudaMemcpyDeviceToDevice));
		if(z_hat != NULL)
			HANDLE_ERROR(cudaMemcpy(z_hat, rhs.z_hat, bytes, cudaMemcpyDeviceToDevice));
	}

	return *this;

}

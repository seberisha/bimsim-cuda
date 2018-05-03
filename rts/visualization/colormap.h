#ifndef RTS_COLORMAP_H
#define RTS_COLORMAP_H

#include <string>
#include <stdlib.h>
#include "rts/cuda/error.h"


#define BREWER_CTRL_PTS 11

void qt_buffer2image(unsigned char* buffer, std::string filename, unsigned int x_size, unsigned int y_size);

static float  BREWERCP[BREWER_CTRL_PTS*4] = {0.192157f, 0.211765f, 0.584314f, 1.0f,
                                      0.270588f, 0.458824f, 0.705882f, 1.0f,
                                      0.454902f, 0.678431f, 0.819608f, 1.0f,
                                      0.670588f, 0.85098f, 0.913725f, 1.0f,
                                      0.878431f, 0.952941f, 0.972549f, 1.0f,
                                      1.0f, 1.0f, 0.74902f, 1.0f,
                                      0.996078f, 0.878431f, 0.564706f, 1.0f,
                                      0.992157f, 0.682353f, 0.380392f, 1.0f,
                                      0.956863f, 0.427451f, 0.262745f, 1.0f,
                                      0.843137f, 0.188235f, 0.152941f, 1.0f,
                                      0.647059f, 0.0f, 0.14902f, 1.0f};


#ifdef __CUDACC__
texture<float4, cudaTextureType1D> cudaTexBrewer;
static cudaArray* gpuBrewer;
#endif



namespace rts{

enum colormapType {cmBrewer, cmGrayscale};

static void buffer2image(unsigned char* buffer, std::string filename, unsigned int x_size, unsigned int y_size)
{
    qt_buffer2image(buffer, filename, x_size, y_size);
}

#ifdef __CUDACC__
static void initBrewer()
{
	//initialize the Brewer colormap

	//allocate CPU space
	float4 cpuColorMap[BREWER_CTRL_PTS];

	//define control rtsPoints
	cpuColorMap[0] = make_float4(0.192157f, 0.211765f, 0.584314f, 1.0f);
	cpuColorMap[1] = make_float4(0.270588f, 0.458824f, 0.705882f, 1.0f);
	cpuColorMap[2] = make_float4(0.454902f, 0.678431f, 0.819608f, 1.0f);
	cpuColorMap[3] = make_float4(0.670588f, 0.85098f, 0.913725f, 1.0f);
	cpuColorMap[4] = make_float4(0.878431f, 0.952941f, 0.972549f, 1.0f);
	cpuColorMap[5] = make_float4(1.0f, 1.0f, 0.74902f, 1.0f);
	cpuColorMap[6] = make_float4(0.996078f, 0.878431f, 0.564706f, 1.0f);
	cpuColorMap[7] = make_float4(0.992157f, 0.682353f, 0.380392f, 1.0f);
	cpuColorMap[8] = make_float4(0.956863f, 0.427451f, 0.262745f, 1.0f);
	cpuColorMap[9] = make_float4(0.843137f, 0.188235f, 0.152941f, 1.0f);
	cpuColorMap[10] = make_float4(0.647059f, 0.0f, 0.14902f, 1.0f);


	int width = BREWER_CTRL_PTS;
	int height = 0;


	// allocate array and copy colormap data
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);

	HANDLE_ERROR(cudaMallocArray(&gpuBrewer, &channelDesc, width, height));

	HANDLE_ERROR(cudaMemcpyToArray(gpuBrewer, 0, 0, cpuColorMap, sizeof(float4)*width, cudaMemcpyHostToDevice));

	// set texture parameters
    cudaTexBrewer.addressMode[0] = cudaAddressModeClamp;
	//texBrewer.addressMode[1] = cudaAddressModeClamp;
    cudaTexBrewer.filterMode = cudaFilterModeLinear;
    cudaTexBrewer.normalized = true;  // access with normalized texture coordinates

	// Bind the array to the texture
    HANDLE_ERROR(cudaBindTextureToArray( cudaTexBrewer, gpuBrewer, channelDesc));

}

static void destroyBrewer()
{
    HANDLE_ERROR(cudaFreeArray(gpuBrewer));

}

template<class T>
__global__ static void applyBrewer(T* gpuSource, unsigned char* gpuDest, unsigned int N, T minVal = 0, T maxVal = 1)
{

	int i = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;

	//compute the normalized value on [minVal maxVal]
	float a = (gpuSource[i] - minVal) / (maxVal - minVal);

	//lookup the color
	float shift = 1.0/(2*BREWER_CTRL_PTS);
	float4 color = tex1D(cudaTexBrewer, a+shift);
	//float4 color = tex1D(cudaTexBrewer, a);

	gpuDest[i * 3 + 0] = 255 * color.x;
	gpuDest[i * 3 + 1] = 255 * color.y;
	gpuDest[i * 3 + 2] = 255 * color.z;
}

template<class T>
__global__ static void applyGrayscale(T* gpuSource, unsigned char* gpuDest, unsigned int N, T minVal = 0, T maxVal = 1)
{
    int i = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;

	//compute the normalized value on [minVal maxVal]
	float a = (gpuSource[i] - minVal) / (maxVal - minVal);

	//threshold
	if(a > 1.0)
        a = 1.0;
    if(a < 0.0)
        a = 0.0;

	gpuDest[i * 3 + 0] = 255 * a;
	gpuDest[i * 3 + 1] = 255 * a;
	gpuDest[i * 3 + 2] = 255 * a;
}

template<class T>
static void gpu2gpu(T* gpuSource, unsigned char* gpuDest, unsigned int nVals, T minVal = 0, T maxVal = 1, colormapType cm = cmGrayscale, int blockDim = 128)
{
	//This function converts a scalar field on the GPU to a color image on the GPU
	int gridX = (nVals + blockDim - 1)/blockDim;
	int gridY = 1;
    if(gridX > 65535)
    {
        gridY = (gridX + 65535 - 1) / 65535;
        gridX = 65535;
    }
    dim3 dimGrid(gridX, gridY);
	//int gridDim = (nVals + blockDim - 1)/blockDim;
	if(cm == cmGrayscale)
		applyGrayscale<<<dimGrid, blockDim>>>(gpuSource, gpuDest, nVals, minVal, maxVal);
	else if(cm == cmBrewer)
	{
		initBrewer();
		applyBrewer<<<dimGrid, blockDim>>>(gpuSource, gpuDest, nVals, minVal, maxVal);
		//HANDLE_ERROR(cudaMemset(gpuDest, 0, sizeof(unsigned char) * nVals * 3));
		destroyBrewer();
	}

}

template<class T>
static void gpu2cpu(T* gpuSource, unsigned char* cpuDest, unsigned int nVals, T minVal, T maxVal, colormapType cm = cmGrayscale)
{
    //this function converts a scalar field on the GPU to a color image on the CPU

    //first create the color image on the GPU

    //allocate GPU memory for the color image
    unsigned char* gpuDest;
    HANDLE_ERROR(cudaMalloc( (void**)&gpuDest, sizeof(unsigned char) * nVals * 3 ));

	//HANDLE_ERROR(cudaMemset(gpuSource, 0, sizeof(T) * nVals));

    //create the image on the gpu
    gpu2gpu(gpuSource, gpuDest, nVals, minVal, maxVal, cm);

	//HANDLE_ERROR(cudaMemset(gpuDest, 0, sizeof(unsigned char) * nVals * 3));

    //copy the image from the GPU to the CPU
    HANDLE_ERROR(cudaMemcpy(cpuDest, gpuDest, sizeof(unsigned char) * nVals * 3, cudaMemcpyDeviceToHost));

	HANDLE_ERROR(cudaFree( gpuDest ));

}

template<typename T>
static void gpu2image(T* gpuSource, std::string fileDest, unsigned int x_size, unsigned int y_size, T valMin, T valMax, colormapType cm = cmGrayscale)
{
	//allocate a color buffer
	unsigned char* cpuBuffer = NULL;
	cpuBuffer = (unsigned char*) malloc(sizeof(unsigned char) * 3 * x_size * y_size);

	//do the mapping
	gpu2cpu<T>(gpuSource, cpuBuffer, x_size * y_size, valMin, valMax, cm);

	//copy the buffer to an image
	buffer2image(cpuBuffer, fileDest, x_size, y_size);

	free(cpuBuffer);
}

#endif

template<class T>
static void cpuApplyBrewer(T* cpuSource, unsigned char* cpuDest, unsigned int N, T minVal = 0, T maxVal = 1)
{
    for(int i=0; i<N; i++)
    {
        //compute the normalized value on [minVal maxVal]
        T v = cpuSource[i];
        float a = (cpuSource[i] - minVal) / (maxVal - minVal);
        if(a < 0) a = 0;
        if(a > 1) a = 1;

        float c = a * (float)(BREWER_CTRL_PTS-1);
        int ptLow = (int)c;
        float m = c - (float)ptLow;
        //std::cout<<m<<std::endl;

        float r, g, b;
        if(ptLow == BREWER_CTRL_PTS - 1)
        {
            r = BREWERCP[ptLow * 4 + 0];
            g = BREWERCP[ptLow * 4 + 1];
            b = BREWERCP[ptLow * 4 + 2];
        }
        else
        {
            r = BREWERCP[ptLow * 4 + 0] * (1.0-m) + BREWERCP[ (ptLow+1) * 4 + 0] * m;
            g = BREWERCP[ptLow * 4 + 1] * (1.0-m) + BREWERCP[ (ptLow+1) * 4 + 1] * m;
            b = BREWERCP[ptLow * 4 + 2] * (1.0-m) + BREWERCP[ (ptLow+1) * 4 + 2] * m;
        }


        cpuDest[i * 3 + 0] = 255 * r;
        cpuDest[i * 3 + 1] = 255 * g;
        cpuDest[i * 3 + 2] = 255 * b;

    }
}

template<class T>
static void cpu2cpu(T* cpuSource, unsigned char* cpuDest, unsigned int nVals, T valMin, T valMax, colormapType cm = cmGrayscale)
{

    if(cm == cmBrewer)
        cpuApplyBrewer(cpuSource, cpuDest, nVals, valMin, valMax);
    else if(cm == cmGrayscale)
    {
        int i;
        float a;
        float range = valMax - valMin;
        for(i = 0; i<nVals; i++)
        {
            //normalize to the range [valMin valMax]
            a = (cpuSource[i] - valMin) / range;

            if(a < 0) a = 0.0;
            if(a > 1) a = 1.0;

            cpuDest[i * 3 + 0] = 255 * a;
            cpuDest[i * 3 + 1] = 255 * a;
            cpuDest[i * 3 + 2] = 255 * a;
        }
    }
}

template<class T>
static void cpu2cpu(T* cpuSource, unsigned char* cpuDest, unsigned int nVals, colormapType cm = cmGrayscale, bool positive = false)
{
    //computes the max and min range automatically

    //find the largest magnitude value
    T maxVal = fabs(cpuSource[0]);
    for(int i=0; i<nVals; i++)
	{
        if(fabs(cpuSource[i]) > maxVal)
            maxVal = fabs(cpuSource[i]);
	}

    if(positive)
        cpu2cpu(cpuSource, cpuDest, nVals, (T)0.0, maxVal, cm);
    else
        cpu2cpu(cpuSource, cpuDest, nVals, -maxVal, maxVal, cm);

}



template<typename T>
static void cpu2image(T* cpuSource, std::string fileDest, unsigned int x_size, unsigned int y_size, T valMin, T valMax, colormapType cm = cmGrayscale)
{
    //allocate a color buffer
	unsigned char* cpuBuffer = (unsigned char*) malloc(sizeof(unsigned char) * 3 * x_size * y_size);

	//do the mapping
	cpu2cpu<T>(cpuSource, cpuBuffer, x_size * y_size, valMin, valMax, cm);

	//copy the buffer to an image
	buffer2image(cpuBuffer, fileDest, x_size, y_size);

	free(cpuBuffer);

}

template<typename T>
static void cpu2image(T* cpuSource, std::string fileDest, unsigned int x_size, unsigned int y_size, colormapType cm = cmGrayscale, bool positive = false)
{
    //allocate a color buffer
	unsigned char* cpuBuffer = (unsigned char*) malloc(sizeof(unsigned char) * 3 * x_size * y_size);

	//do the mapping
	cpu2cpu<T>(cpuSource, cpuBuffer, x_size * y_size, cm, positive);

	//copy the buffer to an image
	buffer2image(cpuBuffer, fileDest, x_size, y_size);

	free(cpuBuffer);

}

}	//end namespace colormap and rts

#endif

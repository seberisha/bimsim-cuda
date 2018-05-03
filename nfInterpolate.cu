#include "nearfield.h"
#include "rts/math/spherical_bessel.h"
#include "rts/math/legendre.h"
#include <stdlib.h>
#include "rts/cuda/error.h"
#include "rts/cuda/timer.h"

texture<float2, cudaTextureType2D> texU_ps1;
texture<float2, cudaTextureType2D> texUf_ps1;


//Incident field for a single plane wave
__global__ void gpuInterpolateU(bsComplex* U, ptype deg, int uR, int vR)
{


    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;

    int uc = uR/2.0;
    int vc = vR/2.0;
    float newx = uc + (iu-uc)*cos(deg) - (iv-vc)*sin(deg) ;
    float newy = vc + (iu-uc)*sin(deg) + (iv-vc)*cos(deg);

//    int newx = uR/2 + ((float)iu-uR/2)*cos(deg) + ((float)iv-vR/2)*sin(deg) ;
//    int newy = vR/2 - ((float)iu-uR/2)*sin(deg) + ((float)iv-vR/2)*cos(deg) ;

    float2 Up = tex2D(texU_ps1, newx + 0.5f, newy + 0.5f);
    U[i] = bsComplex(Up.x, Up.y);
}

//Incident field for a single plane wave
__global__ void gpuInterpolateUf(bsComplex* Uf, ptype deg, int uR, int vR)
{


    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;

    float uc = uR/2.0;
    float vc = vR/2.0;
    float newx = uc + (iu-uc)*cos(deg) - (iv-vc)*sin(deg) ;
    float newy = vc + (iu-uc)*sin(deg) + (iv-vc)*cos(deg);

//    int newx = uR/2 + ((float)iu-uR/2)*cos(deg) + ((float)iv-vR/2)*sin(deg) ;
//    int newy = vR/2 - ((float)iu-uR/2)*sin(deg) + ((float)iv-vR/2)*cos(deg) ;

    float2 Ufp = tex2D(texUf_ps1, newx + 0.5f, newy + 0.5f);

    Uf[i] = bsComplex(Ufp.x, Ufp.y);

}

void nearfieldStruct::saveFieldsToTexture()
{

    //allocate space for the U_ps1 and Uf_ps1 textures
    //allocate the cuda array
    cudaArray* arrayU_ps1;
    cudaArray* arrayUf_ps1;
    cudaChannelFormatDesc channelDescU_ps1 =
            cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
    cudaChannelFormatDesc channelDescUf_ps1 =
            cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);

    HANDLE_ERROR(cudaMallocArray(&arrayU_ps1, &channelDescU_ps1, U.R[0], U.R[1]));
    HANDLE_ERROR(cudaMallocArray(&arrayUf_ps1, &channelDescUf_ps1, Uf.R[0], Uf.R[1]));

    texU_ps1.addressMode[0] = cudaAddressModeBorder;
    texU_ps1.addressMode[1] = cudaAddressModeBorder;
    texU_ps1.filterMode     = cudaFilterModeLinear;
    texU_ps1.normalized     = false;

    texUf_ps1.addressMode[0] = cudaAddressModeBorder;
    texUf_ps1.addressMode[1] = cudaAddressModeBorder;
    texUf_ps1.filterMode     = cudaFilterModeLinear;
    texUf_ps1.normalized     = false;
    HANDLE_ERROR(cudaBindTextureToArray(texU_ps1, arrayU_ps1, channelDescU_ps1));
    HANDLE_ERROR(cudaBindTextureToArray(texUf_ps1, arrayUf_ps1, channelDescUf_ps1));

    //copy U and Uf to texture
    HANDLE_ERROR( cudaMemcpy2DToArray(arrayU_ps1, 0, 0, U.x_hat, U.R[0]*sizeof(float2), U.R[0]*sizeof(float2), U.R[1], cudaMemcpyDeviceToDevice));
    HANDLE_ERROR( cudaMemcpy2DToArray(arrayUf_ps1, 0, 0, Uf.x_hat, Uf.R[0]*sizeof(float2), Uf.R[0]*sizeof(float2), U.R[1], cudaMemcpyDeviceToDevice));

    //gpuStartTimer();


    //t_Uf = gpuStopTimer();
}

void nearfieldStruct::Interpolate()
{
       //#define SQRT_BLOCK 8
    //create one thread for each pixel of the field slice
    dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    dim3 dimGrid((Uf.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (Uf.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);
    //cout<<"\t\t\t rotationAngle = "<<rotationAngle<<endl;
    gpuInterpolateUf<<<dimGrid, dimBlock>>>(Uf.x_hat, rotationAngle, Uf.R[0], Uf.R[1]);
    gpuInterpolateU<<<dimGrid, dimBlock>>>(U.x_hat, rotationAngle, U.R[0], U.R[1]);

}

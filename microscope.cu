#include "microscope.h"

#include "rts/cuda/error.h"
#include "rts/tools/progressbar.h"
#include "rts/cuda/timer.h"
#include "dataTypes.h"
#include "stim/visualization/colormap.h"

#include <stim/image/image.h>

//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
//#include <thrust/generate.h>
//#include <thrust/reduce.h>
//#include <thrust/functional.h>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <sstream>
//#include <sys/time.h>

//#include <QImage>

texture<float, cudaTextureType2D> texUd_ps1;
texture<float, cudaTextureType2D> texUfd_ps1;

texture<float, cudaTextureType1D> texUd_line;
texture<float, cudaTextureType1D> texUfd_line;

extern microscopeStruct* SCOPE;

/*
 * Compute the total field at the detector by interpolating the 1d contribution profile of point
 *  source samples in radial rings
*/
__global__ void InterpolateAndIntegrateUdFromLine(ptype* Ud, int uR, int vR, float scale)
{
    //get the current coordinate in the plane slice
    float iu = blockIdx.x * blockDim.x + threadIdx.x;
    float iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;


    int uc = uR/2.0;
    int vc = vR/2.0;

    float r = sqrt(pow(iu-uc,2.0f) + pow(iv-vc,2.0f));

    Ud[i] = scale*tex1D(texUd_line, r + 0.5f);
}


/*
 * Compute the incident field at the detector by interpolating the 1d contribution profile of point
 *  source samples in radial rings
*/
__global__ void InterpolateAndIntegrateUfdFromLine(ptype* Ufd, int uR, int vR, float scale)
{
    //get the current coordinate in the plane slice
    float iu = blockIdx.x * blockDim.x + threadIdx.x;
    float iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;


    int uc = uR/2.0;
    int vc = vR/2.0;

    float r = sqrt(pow(iu-uc,2.0f) + pow(iv-vc,2.0f));

    Ufd[i] = scale*tex1D(texUfd_line, r + 0.5f);
}


/*
 * Compute the total field at the detector by interpolating the 2d images (at the detector) of point
 *  source samples in radial rings
*/
__global__ void InterpolateAndIntegrateUd(ptype* Ud, ptype deg, int uR, int vR, float scale)
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
    float newy = vc +  (iu-uc)*sin(deg) + (iv-vc)*cos(deg);

    Ud[i] += scale*tex2D(texUd_ps1, newx + 0.5f, newy + 0.5f);
}

/*
 * Compute the incident field at the detector by interpolating the 2d images (at the detector) of point
 *  source samples in radial rings
*/
__global__ void InterpolateAndIntegrateUfd(ptype* Ufd, ptype deg, int uR, int vR, float scale)
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
    Ufd[i] += scale*tex2D(texUfd_ps1, newx + 0.5f, newy + 0.5f);
}

/*
 * Simulate the contribution of a ring of point sources to the final detector image for the incident field.
 * The values at ring samples (at different focal points in the ring) are interpolated from the full computation
 * of detector image of one point source in the ring.
*/

__global__ void simRing_Ufd(ptype* r_d, int uR, int vR, float scale, int numSamples, int numRings)
{

    //shared memory for saving the interp at each sample of the ring
    extern __shared__ float shared[];

    //each block computes interp for a ring
    int radius = blockIdx.x + blockIdx.y * gridDim.x;

    int threadId = threadIdx.x + threadIdx.y*blockDim.x;

    if(radius >= numRings || threadId >= numSamples) return;

    int uc = uR/2.0;
    int vc = vR/2.0;
    float newx = 0.0 ;
    float newy = 0.0;

    //ptype avg = 0.0;
    ptype theta  = (2.0*PI/numSamples)*(threadId+1);

    newx = uc + radius*cos(theta);

    newy = vc - radius*sin(theta);

    shared[threadId] = tex2D(texUfd_ps1, newx + 0.5f, newy + 0.5f);

    //printf("shared[%d], theta %f, %f  \n",threadId, shared[threadId],theta);


    //printf("Thread %d: iter %d\n", threadId, radius);

    // sum reduction

    __syncthreads();

    for(unsigned int step = numSamples/2; step >= 1; step >>= 1)
    {
        __syncthreads();
        if (threadId < step)
        {
            shared[threadId] += shared[threadId + step];
        }
        __syncthreads();
    }
    //__syncthreads();

    if(threadId == 0){
        r_d[radius] += shared[0];///(numSamples + numRings);

        //printf("shared[%d], numSamples: %f, %d\n",radius, shared[0], numSamples);
        //printf("r_d[%d] %f \n",radius, r_d[radius]);
    }

}


// simulate the contribution of a ring of point sources to the final detector image for U
__global__ void simRing_Ud(ptype* r_d, int uR, int vR, float scale, int numSamples, int numRings)
{

    //shared memory for saving the interp at each sample of the ring
    extern __shared__ float shared[];


    //each block computes interp for a ring
    int radius = blockIdx.x + blockIdx.y * gridDim.x;//this is blockID

    int threadId = threadIdx.x + threadIdx.y*blockDim.x;

    if(radius >= numRings || threadId >= numSamples) return;

    int uc = uR/2.0;
    int vc = vR/2.0;
    float newx = 0.0 ;
    float newy = 0.0;

    //ptype avg = 0.0;
    ptype theta  = (2.0*PI/numSamples)*(threadId+1);

    newx = uc + radius*cos(theta);



    newy = vc - radius*sin(theta);

    shared[threadId] = tex2D(texUd_ps1, newx + 0.5f, newy + 0.5f);

    //printf("shared[%d], theta %f, %f  \n",threadId, shared[threadId],theta);


    //printf("Thread %d: iter %d\n", threadId, radius);

    // sum reduction

    __syncthreads();

    for(unsigned int step = numSamples/2; step >= 1; step >>= 1)
    {
        __syncthreads();
        if (threadId < step)
        {
            shared[threadId] += shared[threadId + step];
        }
        __syncthreads();
    }
    //__syncthreads();

    if(threadId == 0){
        r_d[radius] += shared[0];///(numSamples + numRings);

        // printf("Ud \n");
        // printf("r_d[%d] %f \n",radius, r_d[radius]);
    }





    /*
    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    // if(iu > numRings || iv > numRings) return;


    //    if (iu==0 && iv ==0)
    //        printf("\t\t numSamples, numRings = %d,%d\n",numSamples, numRings);

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;

    //    if (i==0)
    //        printf("\t\t numSamples, numRings = %d,%d\n",numSamples, numRings);


    if(i >= numRings) return;

    int uc = uR/2.0;
    int vc = vR/2.0;
    float newx = 0.0 ;
    float newy = 0.0;

    //numSamples = 2*PI*(i+1);

    ptype avg = 0.0;
    ptype theta  = 2.0*PI/numSamples;


    //avg += tex2D(texUfd_ps1, uc + 0.5f, vc + 0.5f);


    // printf("\t\t *** uc,vc %d, %d***\n",uc,vc);
    for (int p = 0; p < numSamples; ++p)
    {
        newx = uc + i*cos(theta);
        newy = vc - i*sin(theta);
        avg += tex2D(texUd_ps1, newx + 0.5f, newy + 0.5f);
        //        if (i==1)
        //        {
        //            printf("\t\t p, ring, avg, iu, iv, newx, newy, theta %d, %d %f, %d, %d, %f, %f, %f\n", p, i, avg, iu,iv,newx,newy,theta);

        //        }
        theta += 2*PI/numSamples;
    }

    r_d[i] += avg;
    */
}



void microscopeStruct::InterpAtDetector()
{
    //#define SQRT_BLOCK 8
    //create one thread for each pixel of the field slice
    dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    dim3 dimGrid((D->R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (D->R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);
    //cout<<"\t\t\t rotationAngle = "<<rotationAngle<<endl;
    InterpolateAndIntegrateUfd<<<dimGrid, dimBlock>>>(Di->S, Di->rotationAngle, Di->R[0], Di->R[1], scale);
    InterpolateAndIntegrateUd<<<dimGrid, dimBlock>>>(D->S, D->rotationAngle, D->R[0], D->R[1], scale);
}


void microscopeStruct::InterpLineAtDetector()
{
    //#define SQRT_BLOCK 8
    //create one thread for each pixel of the field slice
    dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    dim3 dimGrid((D->R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (D->R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

    InterpolateAndIntegrateUfdFromLine<<<dimGrid, dimBlock>>>(Di->S, Di->R[0], Di->R[1], scale);
    InterpolateAndIntegrateUdFromLine<<<dimGrid, dimBlock>>>(D->S, D->R[0], D->R[1], scale);

}

void microscopeStruct::InterpLineForAllPointSources(int numSamples, int numRings)
{
    //#define SQRT_BLOCK 8
    // each block computes a ring
    // each thread in the block computes a ps sample in the ring
    int sqrtThreads = sqrtf(numSamples)+1;
    dim3 dimBlock(sqrtThreads, sqrtThreads);
    int sqrtBlocks = sqrtf(numRings)+1;
    dim3 dimGrid(sqrtBlocks,sqrtBlocks);

    // dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    // dim3 dimGrid((D->R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (D->R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

    simRing_Ufd<<<dimGrid, dimBlock, sizeof(float)*numSamples>>>(Di->line, Di->R[0], Di->R[1], scale, numSamples, numRings);

    simRing_Ud<<<dimGrid, dimBlock, sizeof(float)*numSamples>>>(D->line, D->R[0], D->R[1], scale, numSamples, numRings);

}


void microscopeStruct::SaveDetectorToTexture()
{

    //allocate space for the Ud_ps1 and Ufd_ps1 textures
    //allocate the cuda array
    cudaArray* arrayUd_ps1;
    cudaArray* arrayUfd_ps1;
    cudaChannelFormatDesc channelDescUd_ps1 =
            cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaChannelFormatDesc channelDescUfd_ps1 =
            cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

    HANDLE_ERROR(cudaMallocArray(&arrayUd_ps1, &channelDescUd_ps1, D->R[0], D->R[1]));
    HANDLE_ERROR(cudaMallocArray(&arrayUfd_ps1, &channelDescUfd_ps1, Di->R[0], Di->R[1]));

    texUd_ps1.addressMode[0] = cudaAddressModeBorder;
    texUd_ps1.addressMode[1] = cudaAddressModeBorder;
    texUd_ps1.filterMode     = cudaFilterModeLinear;
    texUd_ps1.normalized     = false;

    texUfd_ps1.addressMode[0] = cudaAddressModeBorder;
    texUfd_ps1.addressMode[1] = cudaAddressModeBorder;
    texUfd_ps1.filterMode     = cudaFilterModeLinear;
    texUfd_ps1.normalized     = false;
    HANDLE_ERROR(cudaBindTextureToArray(texUd_ps1, arrayUd_ps1, channelDescUd_ps1));
    HANDLE_ERROR(cudaBindTextureToArray(texUfd_ps1, arrayUfd_ps1, channelDescUfd_ps1));

    //copy U and Uf to texture
    HANDLE_ERROR( cudaMemcpy2DToArray(arrayUd_ps1, 0, 0, D_ps1->S, D->R[0]*sizeof(float), D->R[0]*sizeof(float), D->R[1], cudaMemcpyDeviceToDevice));
    HANDLE_ERROR( cudaMemcpy2DToArray(arrayUfd_ps1, 0, 0, Di_ps1->S, Di->R[0]*sizeof(float), Di->R[0]*sizeof(float), Di->R[1], cudaMemcpyDeviceToDevice));
}



void microscopeStruct::SaveLineToTexture(int radiusLength)
{

    printf("\t\t SaveLineToTexture\n");
    //allocate space for the U_p and Uf line textures
    //allocate the cuda array
    cudaArray* arrayUd_line;
    cudaArray* arrayUfd_line;
    cudaChannelFormatDesc channelDescUd_line =
            cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaChannelFormatDesc channelDescUfd_line =
            cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    HANDLE_ERROR(cudaMallocArray(&arrayUd_line, &channelDescUd_line, radiusLength, 1));
    HANDLE_ERROR(cudaMallocArray(&arrayUfd_line, &channelDescUfd_line, radiusLength, 1));

    texUd_line.addressMode[0] = cudaAddressModeBorder;
    texUd_line.addressMode[1] = cudaAddressModeBorder;
    texUd_line.filterMode     = cudaFilterModeLinear;
    texUd_line.normalized     = false;

    texUfd_line.addressMode[0] = cudaAddressModeBorder;
    texUfd_line.addressMode[1] = cudaAddressModeBorder;
    texUfd_line.filterMode     = cudaFilterModeLinear;
    texUfd_line.normalized     = false;
    HANDLE_ERROR(cudaBindTextureToArray(texUd_line, arrayUd_line, channelDescUd_line));
    HANDLE_ERROR(cudaBindTextureToArray(texUfd_line, arrayUfd_line, channelDescUfd_line));

    //copy U line and Uf line to texture
    HANDLE_ERROR( cudaMemcpy2DToArray(arrayUd_line, 0, 0, D->line, radiusLength*sizeof(float), radiusLength*sizeof(float), 1, cudaMemcpyDeviceToDevice));
    HANDLE_ERROR( cudaMemcpy2DToArray(arrayUfd_line, 0, 0, Di->line, radiusLength*sizeof(float), radiusLength*sizeof(float), 1, cudaMemcpyDeviceToDevice));
}


__global__ void backpropagate(bsComplex* U, int uR, int vR, ptype du, ptype dv,ptype lambda, float d)
{

    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;

    // backpropagation -- this is done in radians

    float Kx = 2*PI * uR * du ;			//calculate the width of the momentum space
    float Ky = 2*PI * vR * dv;

//    if (i==0)
//        printf("\t\t Kx, Ky, d: %f, %f %f\n",Kx, Ky,d);

    float kx, kx_sq, ky, ky_sq, k_sq;
    float kz;
    bsComplex shift;
    float dkx = Kx / (uR);
    float dky = Ky / (vR);


    /*
     * spatial frequencies if ciruclar shift is performed after fft
    */

    //float min_kx = -Kx / 2;
    //float min_ky = -Ky / 2;

    //kx = min_kx + iu * dkx;					//calculate the position of the current plane wave
    //ky = min_ky + iv * dky;


    if(iu <= uR / 2)
        kx = (ptype)iu * dkx;
    else
        kx = -(ptype)(uR - 1 - iu) * dkx;

    if(iv <= vR / 2)
        ky = (ptype)iv * dky;
    else
        ky = -(ptype)(vR - 1 - iv) * dky;

    kx_sq = kx * kx;
    ky_sq = ky * ky;
    float k = 2*PI/lambda;
    k_sq = k*k;


    if(kx_sq + ky_sq < k_sq){
        kz = sqrt(k*k - kx * kx - ky * ky);			//estimate kz using the Fresnel approximation
        shift = -exp(bsComplex(0, kz * d));
        U[i] *= shift;
        //printf("\t\t shift: %f\n",shift.abs());	
    }
    else{
        //printf("\t\t***************** i = %d\n",i);
        U[i] = 0;
    }

}





__global__ void bandpass(bsComplex* U, int uR, int vR, ptype du, ptype dv, ptype NAin, ptype NAout, ptype lambda)
{

    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;

    ptype u, v;

    if(iu <= uR / 2)
        u = (ptype)iu * du;
    else
        u = -(ptype)(uR - 1 - iu) * du;

    if(iv <= vR / 2)
        v = (ptype)iv * dv;
    else
        v = -(ptype)(vR - 1 - iv) * dv;

    ptype fmag = sqrt(u*u + v*v);
   // printf("\t\t ",);


    if(fmag < NAin / lambda || fmag > NAout / lambda)
        U[i] = 0;
}

microscopeStruct::microscopeStruct()
{
    scalarSim = true;
    D = NULL;
    Di = NULL;
    D_ps1 = NULL;
    Di_ps1 = NULL;
}

void microscopeStruct::init()
{

    nf.scalarSim = scalarSim;
    //Ud.scalarField = scalarSim;
    //Ufd.scalarField = scalarSim;

    //Ud.init_gpu();
    //Ufd.init_gpu();

    //initialize the near field
    nf.init();

    if (SCOPE->interpolate == 2)
    {
        //allocate space for the detector
        D = new scalarslice(Ud.R[0] / ss, Ud.R[1] / ss, SCOPE->nf.extSource[1]);
        Di = new scalarslice(Ud.R[0] / ss, Ud.R[1] / ss, SCOPE->nf.extSource[1]);
    }else
    {
        //allocate space for the detector
        D = new scalarslice(Ud.R[0] / ss, Ud.R[1] / ss);
        Di = new scalarslice(Ud.R[0] / ss, Ud.R[1] / ss);
    }

    D_ps1 = new scalarslice(Ud.R[0] / ss, Ud.R[1] / ss);
    Di_ps1 = new scalarslice(Ud.R[0] / ss, Ud.R[1] / ss);
    //clear the detector
    clearDetector();

}

void microscopeStruct::destroy()
{
    delete D;
    D = NULL;

    delete Di;
    Di = NULL;

    delete D_ps1;
    D_ps1 = NULL;


    delete Di_ps1;
    Di_ps1 = NULL;

    Ud.kill_gpu();
    Ufd.kill_gpu();

    //destroy the near field
    nf.destroy();

}

void microscopeStruct::applyBackpropagation()
{	
    nf.U.toAngularSpectrum();
    //nf.Uf.toAngularSpectrum();
    dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    dim3 dimGrid((nf.U.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (nf.U.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);
    ptype du = 1.0/(nf.pos.X.len());
    ptype dv = 1.0/(nf.pos.Y.len()); 
    //cout<<"\t\t z_d: "<<z_d<<endl; 
    backpropagate<<<dimGrid, dimBlock>>>(nf.U.x_hat, nf.U.R[0], nf.U.R[1], du, dv, nf.lambda, z_d); 
    //backpropagate<<<dimGrid, dimBlock>>>(nf.Uf.x_hat, nf.U.R[0], nf.U.R[1], du, dv, nf.lambda, z_d); 
    nf.U.fromAngularSpectrum();
    //nf.Uf.fromAngularSpectrum();
}


void microscopeStruct::applyBandpass()
{

    //This function applies the objective bandpass to the near field
    //The near field structure stores the results, in order to save memory

    //first convert the near field to an angular spectrum (FFT)
    nf.U.toAngularSpectrum();

    //fft of the focused/incident field
    nf.Uf.toAngularSpectrum();


    //create one thread for each pixel of the field slice
    dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    dim3 dimGrid((nf.U.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (nf.U.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

    //compute the step size in the frequency domain
    ptype du = 1.0 / (nf.pos.X.len());
    ptype dv = 1.0 / (nf.pos.Y.len());

    //apply the objective band-pass filter
    bandpass<<<dimGrid, dimBlock>>>(nf.U.x_hat, nf.U.R[0], nf.U.R[1], du, dv, objective[0], objective[1], nf.lambda);

    //apply bpf to the focused/incident field
    bandpass<<<dimGrid, dimBlock>>>(nf.Uf.x_hat, nf.Uf.R[0], nf.Uf.R[1], du, dv, objective[0], objective[1], nf.lambda);

    nf.U.fromAngularSpectrum();

    nf.Uf.fromAngularSpectrum();

}

void microscopeStruct::getFarField()
{
    //Compute the Far Field image of the focal plane

    //clear the memory from previous detector fields
    //Ud.kill_gpu();
    //Ufd.kill_gpu();

    //first crop the filtered near-field image of the source and scattered fields
    Ud = nf.U.crop(padding * Ud.R[0], padding * Ud.R[1], Ud.R[0], Ud.R[1]);
    //Ud.Mag().toImage("U-ff.bmp");
    Ufd = nf.Uf.crop(padding * Ufd.R[0], padding * Ufd.R[1], Ufd.R[0], Ufd.R[1]);

    //Ufd.Mag().toImage("Uf-ff.bmp");
}

void microscopeStruct::integrateDetector()
{
    scale = 1;
    Ud.IntegrateAndResample(D, ss, scale);
    Ufd.IntegrateAndResample(Di, ss, scale);
}

void microscopeStruct::saveDetectorForOnePs()
{
    //cout<<"\t\t resampling D_ps1"<<endl;
    Ud.ResampleAndSave(D_ps1, ss);
    //cout<<"\t\t resampling Di_ps1"<<endl;
    Ufd.ResampleAndSave(Di_ps1, ss);
}

void microscopeStruct::clearDetector()
{
    //zero-out the detector
    D->clear();
    Di->clear();
    //clear these only if using interp with images of point sources
    D_ps1->clear();
    Di_ps1->clear();
    if (SCOPE->interpolate == 2)
    {
        Di->clearLine(SCOPE->nf.extSource[1]);
        D->clearLine(SCOPE->nf.extSource[1]);
    }

}

//flag for a vector simulation
void microscopeStruct::setPos(bsPoint pMin, bsPoint pMax, bsVector normal)
{
    pos = rts::quad<ptype, 3>(pMin, pMax, normal);
}

//set slice distance from the center of the sphere
void microscopeStruct::setZdPos(float d)
{
    z_d = d;
}

void microscopeStruct::setRes(int x_res, int y_res, int pad, int supersampling)
{
    padding = pad;
    ss = supersampling;

    Ufd.R[0] = Ud.R[0] = x_res * ss;
    Ufd.R[1] = Ud.R[1] = y_res * ss;
}

void microscopeStruct::setNearfield()
{
    //sets the values for the near field in order to create the specified detector image

    //compute the size of the near-field slice necessary to create the detector image
    nf.pos = pos * (padding * 2 + 1);
    //cout<<"\t\t nf.pos, padding "<<nf.pos.toStr()<<" , "<<padding<<endl;

    //compute the resolution of the near-field slice necessary to create the detector image
    nf.setRes(Ud.R[0] * (padding * 2 + 1), Ud.R[1] * (padding * 2 + 1));

}

__global__ void calc_absorbance(ptype* A, ptype* D, ptype* Di, int N)
{
    //compute the index for this thread
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;

    A[i] = -log10(D[i] / Di[i]);
    //printf("D[i], Di[i], A[i] %f , %f , %f\n", D[i], Di[i], A[i]);

}

scalarslice microscopeStruct::getAbsorbance()
{
    //compute the magnitude of the field at each rtsPoint in the slice

    //create a scalar slice at the same resolution as the field
    scalarslice* A = new scalarslice(D->R[0], D->R[1]);

    //compute the total number of values in the slice
    unsigned int N = D->R[0] * D->R[1];
    int gridDim = (N+BLOCK-1)/BLOCK;

    calc_absorbance<<<gridDim, BLOCK>>>(A->S, D->S, Di->S, N);

    return *A;
}

ptype microscopeStruct::getAbsorbanceSpectrum()
{

    //compute the absorbance spectrum at a particular wavelength

    //thrust::device_ptr<ptype> thD_ptr(D->S);
    //ptype sumD = thrust::reduce(thD_ptr, thD_ptr + (D->R[0]*D->R[1]));
    //thrust::device_ptr<ptype> thDi_ptr(Di->S);
    //ptype sumDi = thrust::reduce(thDi_ptr,thDi_ptr + (Di->R[0]*Di->R[1]));

    //return -log10(sumD/sumDi);
    return 0;

}



__global__ void calc_transmittance(ptype* A, ptype* D, ptype* Di, int N)
{
    //compute the index for this thread
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;


    A[i] = D[i] / Di[i];

}

scalarslice microscopeStruct::getTransmittance()
{
    //compute the magnitude of the field at each rtsPoint in the slice

    //create a scalar slice at the same resolution as the field
    scalarslice* T = new scalarslice(D->R[0], D->R[1]);

    //compute the total number of values in the slice
    unsigned int N = D->R[0] * D->R[1];
    int gridDim = (N+BLOCK-1)/BLOCK;

    calc_transmittance<<<gridDim, BLOCK>>>(T->S, D->S, Di->S, N);

    return *T;
}

scalarslice microscopeStruct::getIntensity()
{
    //create a scalar slice at the same resolution as the field
    scalarslice* I = new scalarslice(D->R[0], D->R[1]);

    HANDLE_ERROR(cudaMemcpy(I->S, D->S, sizeof(ptype) * D->R[0] * D->R[1], cudaMemcpyDeviceToDevice));

    return *I;

}


scalarslice microscopeStruct::getIncidentFieldImage()
{
    //create a scalar slice at the same resolution as the field
    scalarslice* I = new scalarslice(D->R[0], D->R[1]);

    HANDLE_ERROR(cudaMemcpy(I->S, Di->S, sizeof(ptype) * D->R[0] * D->R[1], cudaMemcpyDeviceToDevice));

    return *I;

}



void microscopeStruct::SimulateScattering()
{
    //cout<<"\t\t focal point: "<<nf.focus.toStr()<<endl;
    //cout<<"\t\t amplitude: "<<nf.A<<endl;
    nf.Simulate();

}


void microscopeStruct::saveToTexture()
{
    //save fields to texture
    nf.saveFieldsToTexture();
}


void microscopeStruct::SimulateInterpolation()
{
    //interpolate fields
    nf.Interpolate();
}

void microscopeStruct::SimulateImaging()
{
    applyBandpass();

    //std::cout<<"\t\t d: "<<z_d<<endl<<endl;
    getFarField();
    integrateDetector();
}

void microscopeStruct::Simulate()
{
    //cout<<"\t\t microscopeStruct::Simulate()"<<endl;
    SimulateScattering();
    //applyBackpropagation();
    SimulateImaging();
    //scalarslice I = getAbsorbance();
    //I.toImage("absImg.bmp");
}


void microscopeStruct::SimulateAndSaveTexture()
{
    SimulateScattering();
    saveToTexture();
    SimulateImaging();
}

void microscopeStruct::SimulateWithInterpolation()
{
    SimulateInterpolation();
    SimulateImaging();
}

void microscopeStruct::SimulateExtendedSource()
{

    clearDetector();
    scale = 1.0;


    //for each source in the source list
    int npts = focalPoints.size();
    //std::cout<<"\t\t npts "<<npts<<std::endl;
    float t=0;
    for(int i = 0; i<npts; i++)
    {
        nf.focus = focalPoints[i].f;
        //cout<<"microscopeStruct::SimulateExtendedSource===nf.focus "<<nf.focus[0]<<", "<<nf.focus[1]<<", "<<nf.focus[2]<<endl;
        nf.A = focalPoints[i].A;

        gpuStartTimer();
        Simulate();
        /*
    	ostringstream tempFile;
        tempFile<<"out_Ufd"<<i<<"_"<<0<<".bmp";
        scalarslice I = getIncidentFieldImage();
        I.toImage(tempFile.str());
    	*/

        //if(verbose)
        //  saveDetectorForEachPointSource();
        t += gpuStopTimer();

        rtsProgressBar((double)(i+1)/(double)npts * 100);
    }
    if(verbose)
    {
        cout<<endl;
        cout<<"Total time: "<<t<<"ms"<<endl;
        cout<<"Time per source: "<<t/npts<<"ms"<<endl;

    }

}


void microscopeStruct::saveDetectorForEachPointSource()
{
    string intFile = "tempU.bmp";
    string incFile = "tempUf.bmp";
    string absFile = "tempA.bmp";

    //intensity
    if(intFile != "")
    {
        scalarslice I = getIntensity();


        I.toImage(intFile);
    }

    //incident field image

    //intensity
    if(incFile != "")
    {
        scalarslice I = getIncidentFieldImage();
        //  std::cout<<"\t\t=====> incident file: "<<incFile<< std::endl;


        I.toImage(incFile);
    }


    //std::cout<<"\t\t******=====> writing absorbance file: "<<std::endl;
    //absorbance
    if(absFile != "")
    {
        scalarslice I = getAbsorbance();

        I.toImage(absFile);
    }

    //std::cout<<"\t\t******=====> writing absorbance spectrum file: "<<std::endl;

    //absorbance
    //    if(absSpecFile != "")
    //    {
    //        ptype a = scope->getAbsorbanceSpectrum();

    //        // std::cout<<"\t\t=====> a is: "<<a<<std::endl;



    //        if(is_binary(absSpecFile))
    //        {
    //            ofstream myfile;
    //            myfile.open (absSpecFile.c_str(), ios::app);

    //            if(wavenumber)
    //                myfile << 10000.0f/scope->nf.lambda <<","<<a<<endl;
    //            else
    //                myfile << scope->nf.lambda <<","<<a<<endl;
    //            myfile.close();


    //        }
    //        //  else
    //        //  I.toImage(absSpecFile);
    //    }


    //    //transmittance
    //    if(transFile != "")
    //    {
    //        scalarslice I = scope->getTransmittance();

    //        if(is_binary(transFile))
    //        {
    //            if(wavenumber)
    //                I.toEnvi(transFile, 10000.0f/scope->nf.lambda, append);
    //            else
    //                I.toEnvi(transFile, scope->nf.lambda, append);
    //        }
    //        else
    //            I.toImage(transFile);
    //    }

}





void microscopeStruct::SimulateExtendedSourceInterpAtDetector(int *pointsRings)
{

    cout<<"\t\t ***Interpolating images of point sources at the detector***"<<endl;

    gpuStartTimer();

    clearDetector();

    int numRings = pointsRings[1];
    float t=0;

    scale = 1.0;    //scale intensity at the detector

    //simulate for point source at fp = (0,0,0)
    Simulate();

    int npts = pointsRings[0];

    float spacing = 1;
    float radius = spacing;
    float amplitude = 1;
    int psIdx = 1;

    int totalPS = npts*(numRings-1) + 1;
    //generate point sources in rings and simulate
    for(float j=spacing; j<numRings; j+=spacing)
    {
        // generate set of point sources in a ring at a particular radius
        generatePointSources(npts, radius, amplitude);

        //compute the full field simulation for the first point source
        nf.focus = focalPointsES[0].f;

        nf.A = focalPointsES[0].A;

        // full field simulation for one point source in the ring
        Simulate();
        ++psIdx;

        /*
        ostringstream tempFile;
        tempFile<<"out_Ufd"<<j<<"_"<<0<<".bmp";
        scalarslice I = getIncidentFieldImage();
        I.toImage(tempFile.str());
        */

        saveDetectorForOnePs();

        SaveDetectorToTexture();

        //obtain the rest of point sources in the same ring by interpolating the first point source
        for(int i = 1; i<npts; ++i)
        {
            nf.focus = focalPointsES[i].f;
            //cout<<"\t\t other focal points = "<<nf.focus.toStr()<<endl;
            nf.A = focalPointsES[i].A;
            nf.rotationAngle = focalPointsES[i].rotationAngle;
            D->rotationAngle = focalPointsES[i].rotationAngle;
            Di->rotationAngle = focalPointsES[i].rotationAngle;
            InterpAtDetector();
            ++psIdx;
            //tempFile = tempFile + to_string(i) + ".bmp";

            /*
             * this prints the 'total' focused field -- not each interpolated one
            ostringstream tempFile;
            tempFile<<"out_Ufd"<<j<<"_"<<i<<".bmp";
            scalarslice I = getIncidentFieldImage();
            I.toImage(tempFile.str());
            */
            rtsProgressBar((double)(psIdx)/(double)totalPS * 100);
        }

        radius += spacing;
        //cout<<"\t\t================= "<<j<<endl;
    }

    t += gpuStopTimer();

    if(verbose)
    {
        cout<<endl;
        cout<<"Total time: "<<t<<endl;
        cout<<"Total number of point source samples: "<< totalPS<<endl;
        cout<<"Time per source: "<<t/totalPS<<"ms"<<endl;
    }

}




void microscopeStruct::SimulateExtendedSourceInterpLineAtDetector(int *pointsRings)
{

    cout<<"\t\t ***Interpolating a line at the detector***"<<endl;

    /* Generate point sources in forms of rings around the sphere.
     * Simulate an extended source computation using interpolation.
     */


    clearDetector();

    int numRings = pointsRings[1]-1;

    //cout<<"\t\t numSamples, numRings: "<<pointsRings[0]<<" , "<<numRings<<endl;
    float t=0;

    scale = 1.0;
    nf.A = 1;

    gpuStartTimer();


    //simulate for the point source at [0 0 0]
    Simulate();

    saveDetectorForOnePs();

    //    thrust::device_ptr<ptype> thD_ptr(D->S);
    //    ptype sumD = thrust::reduce(thD_ptr, thD_ptr + (D->R[0]*D->R[1]));

    //save detector image to texture
    SaveDetectorToTexture();
    //compute the contribution of this point source to the radius of the final detector image
    InterpLineForAllPointSources(pointsRings[0], pointsRings[1]);

    float spacing = 1;

    int numPoints = 0;

    numPoints = pointsRings[0];


    scale = 1.0;

    //simulate for the rest of the point source rings using one line interpolation
    for(int j=spacing; j<=numRings; j+=spacing)
    {
        //gausAtR = exp(-(pow(j,2))/(2*pow(sigma,2)))/(2*PI*pow(sigma,2));
        //amplitude = gausAtR/maxGaus;


        //focal point for the first point source in the ring
        nf.focus[0] = j;
        nf.focus[1] = 0;
        nf.focus[2] = 0;
        nf.A = 1;

        //Ud.clear_gpu();
        // Ufd.clear_gpu();
        //compute the full field simulation for the first point source in the ring
        Simulate();

        saveDetectorForOnePs();

        SaveDetectorToTexture();

        InterpLineForAllPointSources(pointsRings[0], pointsRings[1]);
    }

    SaveLineToTexture(pointsRings[1]);
    //D->clear();
    //Di->clear();
    InterpLineAtDetector();

    t += gpuStopTimer();


    if(verbose)
    {
        cout<<endl;
        cout<<"Total time: "<<t<<endl;
        cout<<"Total number of point source samples: "<<numPoints*(numRings)+1<<endl;
        cout<<"Time per source: "<<t/(numPoints*(numRings)+1)<<"ms"<<endl;
    }

}





void microscopeStruct::SimulateExtendedSourceWithInterpolation(float sphereRadius, int *pointsRings)
{
    /* Generate point sources in forms of rings around the sphere.
     * Simulate an extended source computation using interpolation.
     */
    scale = 1;
    float numPoints = pointsRings[0];
    cout<<"\t\t***SimulateExtendedSourceWithInterpolation***"<<endl;
    clearDetector();
    //cout<<"\t\t\t points rings = "<<pointsRings[0]<<" , "<<pointsRings[1]<<endl;

    int numRings = pointsRings[1];
    float t=0;

    //simulate for the point source at [0 0 0]
    gpuStartTimer();
    Simulate();
    int npts=0;

    float spacing = 1;
    float radius = spacing;
    float amplitude = 0.0;

    int totalPS = numPoints*(numRings-1)+1;
    int psIdx = 1;

    amplitude = 1;

    //generate point sources in rings and simulate
    for(int j=spacing; j<numRings; j+=spacing)
    {

        npts = pointsRings[0];
        generatePointSources(npts, radius, amplitude);
        radius += spacing;

        nf.focus = focalPointsES[0].f;
        nf.A = focalPointsES[0].A;
        gpuStartTimer();
        SimulateAndSaveTexture();
        ++psIdx;

        //interpolate for the rest of the point sources
        for(int i = 1; i<npts; ++i)
        {
            nf.focus = focalPointsES[i].f;
            //cout<<"\t\t other focal points = "<<nf.focus.toStr()<<endl;
            nf.A = focalPointsES[i].A;

            //cout<<"\t\t amplitude: "<<nf.A<<endl;
            nf.rotationAngle = focalPointsES[i].rotationAngle;
            gpuStartTimer();
            SimulateWithInterpolation();
            ++psIdx;
            t += gpuStopTimer();
            rtsProgressBar((double)(psIdx)/(double)totalPS * 100);
        }
    }

    if(verbose)
    {
        cout<<endl;
        cout<<"Total time: "<<t<<endl;
        cout<<"Total number of point source samples: "<<totalPS<<endl;
        cout<<"Time per source: "<<t/(totalPS)<<"ms"<<endl;
    }

}

void microscopeStruct::LoadExtendedSource(std::string filename)
{
    //this function loads an image of an extended source and creates a list of corresponding point sources

    stim::image<unsigned char> sourceImage(filename.c_str());

    //get the resolution of the image (the boundary is scaled to match the detector)
    int Rx = sourceImage.width();

    int Ry = sourceImage.height();

    //for each pixel
    int x, y;
    float u, v;

    for(x=0; x<Rx; x++)
        for(y=0; y<Ry; y++)
        {
            //create a new point source
            sourcePoint p;

            //compute the coordinate of the focal point
            u = (ptype)x / (ptype)Rx;
            v = (ptype)y / (ptype)Ry;
            p.f = pos(u, v);
            //p.f[2] = 0;
            //get the amplitude of the focal point
            unsigned char s = sourceImage(x, y, 0);
            if(s != 0)
            {
                p.A = (ptype) s / 255;
                focalPoints.push_back(p);
           	//cout<<"\t\t"<<p.f.toStr()<<endl;
	    }
        }

    if (verbose)
        cout<<"Number of point sources: "<<focalPoints.size()<<endl;
}

void microscopeStruct::generatePointSources(int numPoints, float radius, float amplitude)
{
    //create a new point source
    sourcePointES p;
    ptype u, v;
    focalPointsES.clear();

    float theta = 2*PI/(numPoints);
    for(int i = 0; i < numPoints; ++i)
    {
        p.rotationAngle = theta;

        //most of these are not used -- for interpolation only need to know the coordinates
        //on the radius, we need the rotation angle though
        //pol2cart
        u = (cos(theta)*radius);// + halfFOV)/(float)(2*halfFOV);//fovSize;
        v = (sin(theta)*radius );//+ halfFOV)/(float)(2*halfFOV);//fovSize;

        p.f[0] = u;
        p.f[1] = v;
        p.f[2] = 0;
        p.A = amplitude;
        focalPointsES.push_back(p);
        theta += 2*PI/(numPoints);

    }
}


void microscopeStruct::SimulateExtendedGausssianSourceNoInterp(int *pointsRings)
{

    cout<<"\t\t ***SimulateExtendedGausssianSourceNoInterp***"<<endl;

    /* Generate point sources in forms of rings around the sphere.
     * Simulate an extended source computation without interpolation.
     * This function performs the computations for point sources in form of rings
     * with Gaussian amplitudes.
     */

    float sigma = pointsRings[0];
    clearDetector();

    int numRings = pointsRings[1];
    float t=0;

    //simulate for the point source at [0 0 a]
    gpuStartTimer();
    Simulate();
    t += gpuStopTimer();

    float spacing = 0.5;
    float radius = spacing;
    float amplitude = 0.0;

    float maxGaus = exp(-(pow(0,2))/(2*pow(sigma,2)))/(2*PI*pow(sigma,2));
    cout<<"\t\t maxGaus: "<<maxGaus<<endl;
    float minGaus =0.0;

    if (numRings==1)
        minGaus = exp(-(pow(0,2))/(2*pow(sigma,2)))/(2*PI*pow(sigma,2));
    else
        minGaus = exp(-(pow(numRings-1,2))/(2*pow(sigma,2)))/(2*PI*pow(sigma,2));


    float oldRange = maxGaus - minGaus;

    float newRange = 1.0;
    float gausAtR = 0.0;


    //generate point sources in rings and simulate
    for(int j=spacing; j<numRings; j+=spacing)
    {
        gausAtR = exp(-(pow(j,2))/(2*pow(sigma,2)))/(2*PI*pow(sigma,2));
        amplitude =  (((gausAtR - minGaus) * newRange) / oldRange);

        //if (amplitude > 1e-7)
        //{
        generatePointSourcesNoInterp(spacing, radius, amplitude);

        radius += spacing;

        //compute the full field simulation for the first point source

        gpuStartTimer();

        //for each source in the source list
        int npts = focalPoints.size();
        float t=0;
        for(int i = 0; i<npts; i++)
        {
            nf.focus = focalPoints[i].f;
            nf.A = focalPoints[i].A;

            gpuStartTimer();
            Simulate();

            t += gpuStopTimer();

        }
        if(verbose)
        {
            cout<<endl;
            cout<<"Time per source: "<<t/npts<<"ms"<<endl;
        }


        //}
    }

}


void microscopeStruct::generatePointSourcesNoInterp(float  spacing, int radius, float amplitude)
{
    /*generate point sources in rings for angles between 0 and 2pi for the
      no interpolation version -- gaussian source
    */

    //create a new point source
    sourcePoint p;
    ptype u, v;
    focalPoints.clear();
    //    cout<<"\t\t\t sphereRadius = "<<sphereRadius<<endl;
    // cout<<"\t\t\t spacing = "<<spacing<<endl;


    float numPoints = std::round((2*PI*radius-spacing)/spacing + 1);
    int maxRange = std::round(2*PI*radius - spacing);
    float s = 0;



    if(numPoints==1)
    {
        u = 0;
        v = 0;
        p.f[0] = 0;
        p.f[1] = 0;
        p.f[2] = 0;
        p.A = 1;
        focalPoints.push_back(p);

    }else
    {

        ptype a = 0, b = maxRange;

        ptype step  = (b-a)/(numPoints-1);

        for(int i = 0; i < numPoints; ++i)
        {

            a = s/radius;
            u = cos(a)*radius;//Ud.R[0];
            v = sin(a)*radius;//Ud.R[1];

            p.f[0] = u;
            p.f[1] = v;
            p.f[2] = 0;

            p.A = amplitude;
            focalPoints.push_back(p);
            s = s + step;
        }
    }
}


std::string microscopeStruct::toStr()
{
    stringstream ss;
    ss<<nf.toStr();

    ss<<"----------Optics--------------"<<endl<<endl;
    ss<<"Objective NA: "<<objective[0]<<" to "<<objective[1]<<endl;
    return ss.str();


}

#include "nearfield.h"
#include "rts/math/spherical_bessel.h"
#include "rts/math/legendre.h"
#include <stdlib.h>
#include "rts/cuda/error.h"
#include "rts/cuda/timer.h"

texture<float2, cudaTextureType2D> texUsp;
texture<float2, cudaTextureType2D> texUip;

__global__ void gpuScalarUpLut(bsComplex* Us, bsVector* k, int nk, ptype kmag, ptype a, ptype dmin, ptype dmax, bsPoint f, bsPoint ps, ptype A, bsRect ABCD, int uR, int vR, int dR, int aR, int thetaR)
{
    /*This function uses Monte-Carlo integration to sample a texture-based LUT describing the scattered field
        produced by a plane wave through a sphere.  The MC sampling is used to approximate a focused field.

        Us  =   final scattered field
        k   =   list of incoming plane waves (Monte-Carlo samples)
        nk  =   number of incoming MC samples
        kmag=   magnitude of the incoming field 2pi/lambda
        dmin=   minimum distance of the Usp texture
        dmax=   maximum distance of the Usp texture
        f   =   position of the focus
        ps  =   position of the sphere
        A   =   total amplitude of the incident field arriving at the focal spot
        ABCD=   rectangle representing the field slice
        uR  =   resolution of the field slice in the u direction
        vR  =   resolution of the field slice in the v direction
        dR  =   resolution of the Usp texture in the d direction
        thetaR= resolution of the Usp texture in the theta direction
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
    bsVector r = p - ps;
    bsVector distance = p - f;
    float normDistance = distance.len();

    ptype d = r.len();
    float di = ( (d - max(a, dmin))/(dmax - max(a, dmin)) ) * (dR - 1);
    float ai = ( (d - dmin)/(a - dmin)) * (aR - 1);

    bsComplex sumUs(0, 0);
    //for each plane wave in the wave list
    for(int iw = 0; iw < nk; iw++)
    {
        //normalize the direction vectors and find their inner product
        r = r.norm();
        ptype cos_theta = k[iw].dot(r);
        if(cos_theta < -1)
            cos_theta = -1;
        if(cos_theta > 1)
            cos_theta = 1;
        float thetai = ( acos(cos_theta) / PI ) * (thetaR - 1);

        //compute the phase factor for spheres that are not at the origin
        bsVector c = ps - f;
        bsComplex phase = exp(bsComplex(0, kmag * k[iw].dot(c)));
       // bsVector sliceVec(0,0,3.25);
        //bsVector z = sliceVec;
        //ptype = sqrt(pow(kmag,2) - pow(p[0],2) - pow([1],2));

       // bsComplex green = exp(bsComplex(0, kmag * k[iw].dot(z)));
//        if(iu==0 && iv==0)
//            printf("\t\t green %f,%f, %f\n", green.real(), green.imag(), green.abs());
        //compute the internal field if we are inside a sphere
        if(d < a)
        {
            //printf("\t\t====================== internal field \n");
            float2 Uip = tex2D(texUip, ai + 0.5f, thetai + 0.5f);
            sumUs += (1.0f/nk) * A * phase * bsComplex(Uip.x, Uip.y);
        }
        //otherwise compute the scattered field
        else
        {
            float2 Usp = tex2D(texUsp, di + 0.5f, thetai + 0.5f);
            sumUs += (1.0f/nk) * A *phase * bsComplex(Usp.x, Usp.y);
        }

    }

    //Us[i] += sumUs/normDistance;
    Us[i] += sumUs;
}



__global__ void opt_gpuScalarUpLut(bsComplex* Us, bsVector* k, int nk, ptype kmag, ptype a, ptype dmin, ptype dmax, bsPoint f, bsPoint ps, ptype A, bsRect ABCD, int uR, int vR, int dR, int aR, int thetaR)
{
    /*This function uses Monte-Carlo integration to sample a texture-based LUT describing the scattered field
        produced by a plane wave through a sphere.  The MC sampling is used to approximate a focused field.

        Us  =   final scattered field
        k   =   list of incoming plane waves (Monte-Carlo samples)
        nk  =   number of incoming MC samples
        kmag=   magnitude of the incoming field 2pi/lambda
        dmin=   minimum distance of the Usp texture
        dmax=   maximum distance of the Usp texture
        f   =   position of the focus
        ps  =   position of the sphere
        A   =   total amplitude of the incident field arriving at the focal spot
        ABCD=   rectangle representing the field slice
        uR  =   resolution of the field slice in the u direction
        vR  =   resolution of the field slice in the v direction
        dR  =   resolution of the Usp texture in the d direction
        thetaR= resolution of the Usp texture in the theta direction
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
    bsVector r = p - ps;
    ptype d = r.len();
    float di = ( (d - max(a, dmin))/(dmax - max(a, dmin)) ) * (dR - 1);
    float ai = ( (d - dmin)/(a - dmin)) * (aR - 1);

    //normalize the direction vectors and find their inner product
    r = r.norm();


    bsComplex sumUs(0, 0);
    //compute the phase factor for spheres that are not at the origin
    bsVector c = ps - f;


    if(d < a)
    {
        for(int iw = 0; iw < nk; iw++)
        {
            ptype cos_theta = k[iw].dot(r);
            if(cos_theta < -1)
                cos_theta = -1;
            if(cos_theta > 1)
                cos_theta = 1;
            float thetai = ( acos(cos_theta) / PI ) * (thetaR - 1);

            bsComplex phase = exp(bsComplex(0, kmag * k[iw].dot(c)));


            float2 Uip = tex2D(texUip, ai + 0.5f, thetai + 0.5f);
            sumUs += (1.0f/nk) * A * phase * bsComplex(Uip.x, Uip.y);
        }

    }else
    {   for(int iw = 0; iw < nk; iw++)
        {
            ptype cos_theta = k[iw].dot(r);
            if(cos_theta < -1)
                cos_theta = -1;
            if(cos_theta > 1)
                cos_theta = 1;
            float thetai = ( acos(cos_theta) / PI ) * (thetaR - 1);

            bsComplex phase = exp(bsComplex(0, kmag * k[iw].dot(c)));


            float2 Usp = tex2D(texUsp, di + 0.5f, thetai + 0.5f);
            sumUs += (1.0f/nk) * A * phase * bsComplex(Usp.x, Usp.y);
        }
    }

    Us[i] += sumUs;
}



__global__ void optPc_gpuScalarUpLut(ptype* d, bsComplex* Us, int nk,ptype a,
                                     int uR, int vR, ptype *ai, ptype *di, ptype *thetai,
                                     bsComplex *scaledPhase)
{

    /*This function uses Monte-Carlo integration to sample a texture-based LUT describing the scattered field
        produced by a plane wave through a sphere.  The MC sampling is used to approximate a focused field.

        Us  =   final scattered field
        k   =   list of incoming plane waves (Monte-Carlo samples)
        nk  =   number of incoming MC samples
        kmag=   magnitude of the incoming field 2pi/lambda
        dmin=   minimum distance of the Usp texture
        dmax=   maximum distance of the Usp texture
        f   =   position of the focus
        ps  =   position of the sphere
        A   =   total amplitude of the incident field arriving at the focal spot
        ABCD=   rectangle representing the field slice
        uR  =   resolution of the field slice in the u direction
        vR  =   resolution of the field slice in the v direction
        dR  =   resolution of the Usp texture in the d direction
        thetaR= resolution of the Usp texture in the theta direction
    */



    /****use shared memory for scaled phase***/
    extern __shared__ bsComplex sh_scaledPhase[];

    int threadId = threadIdx.y*blockDim.x + threadIdx.x;

    int blockSize = blockDim.x*blockDim.y;

    if (blockSize < nk)
    {
        int endIdx = (int) ceilf((float)nk/blockSize);

        for (int i = 0;i < endIdx; ++i)
            sh_scaledPhase[threadId+i*blockSize] = scaledPhase[threadId+i*blockSize];

    }else
    {

        if (threadId < nk)
            sh_scaledPhase[threadId] = scaledPhase[threadId];

    }

    __syncthreads();



    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;
    //int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    //  int threadId = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;



    bsComplex sumUs(0, 0);


    //printf("\t\t\t nk is %d\n", nk);

    for(int iw = 0; iw < nk; iw++)
    {
        if(d[i] < a)
        {

            float2 Uip = tex2D(texUip, ai[i] + 0.5f,  thetai[iw*uR*vR + i] + 0.5f);
            sumUs += sh_scaledPhase[iw] * bsComplex(Uip.x, Uip.y);
        }else
        {
            float2 Usp = tex2D(texUsp, di[i] + 0.5f, thetai[iw*uR*vR + i] + 0.5f);
            sumUs += sh_scaledPhase[iw] * bsComplex(Usp.x, Usp.y);
        }

    }



    Us[i] += sumUs;
}




__global__ void precomputeLUT(int dR, int aR, ptype dmax,ptype dmin, bsRect ABCD, int nk,
                              ptype *thetai, bsVector* k,
                              ptype *ai, ptype *di, ptype *d, bsPoint ps, int uR, int vR,
                              int thetaR, ptype a, bsComplex *scaledPhase, ptype kmag, bsVector c,
                              ptype A)
{


    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;

    if (i < nk){

        scaledPhase[i] = (1.0f/nk) * A * exp(bsComplex(0, kmag * k[i].dot(c)));
        //  printf("\t\t scaledPhase[%d]: %f\n", i, scaledPhase[i].abs());
    }

    //compute the parameters for u and v
    ptype u = (ptype)iu / (uR);
    ptype v = (ptype)iv / (vR);

    //get the rtsPoint in world space and then the r vector
    bsPoint p = ABCD(u, v);
    bsVector r = p - ps;
    d[i] = r.len();

    di[i] = ( (d[i] - max(a, dmin))/(dmax - max(a, dmin)) ) * (dR - 1);
    ai[i] = ( (d[i] - dmin)/(a - dmin)) * (aR - 1);



    //normalize the direction vectors and find their inner product
    r = r.norm();
    for(int iw = 0; iw < nk; iw++)
    {
        ptype cos_theta = k[iw].dot(r);
        if(cos_theta < -1)
            cos_theta = -1;

        if(cos_theta > 1)
            cos_theta = 1;
        thetai[iw*uR*vR + i] = ( acos(cos_theta) / PI ) * (thetaR - 1);

    }

}


void nearfieldStruct::scalarUpLut()
{
    //get the number of spheres
    int nSpheres = sVector.size();

    //if there are no spheres, nothing to do here
    if(nSpheres == 0)
        return;

    //time the calculation of the focused field
    gpuStartTimer();

    //clear the scattered field
    U.clear_gpu();

    //create one thread for each pixel of the field slice
    dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    dim3 dimGrid((U.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (U.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

    //copy Monte-Carlo samples to the GPU and determine the incident amplitude (plane-wave specific stuff)
    bsVector* gpuk;
    int nWaves;
    ptype subA;
    if(planeWave)
    {
        nWaves = 1;
        HANDLE_ERROR(cudaMalloc( (void**)&gpuk, sizeof(bsVector) ) );
        HANDLE_ERROR(cudaMemcpy( gpuk, &k, sizeof(bsVector), cudaMemcpyHostToDevice));
        subA = A;
    }
    else
    {
        nWaves = inWaves.size();
        HANDLE_ERROR(cudaMalloc( (void**)&gpuk, sizeof(bsVector) * nWaves ) );
        HANDLE_ERROR(cudaMemcpy( gpuk, &inWaves[0], sizeof(bsVector) * nWaves, cudaMemcpyHostToDevice));
        //compute the amplitude that makes it through the condenser
        subA = 2 * PI * A * ( (1 - cos(asin(condenser[1]))) - (1 - cos(asin(condenser[0]))) );
    }

    //for each sphere
    for(int s = 0; s<nSpheres; s++)
    {
        //get the current sphere
        //sphere S = sVector[s];

        //allocate space for the Usp and Uip textures
        //allocate the cuda array
        cudaArray* arrayUsp;
        cudaArray* arrayUip;
        cudaChannelFormatDesc channelDescUsp =
                cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
        cudaChannelFormatDesc channelDescUip =
                cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
        int dR = sVector[s].Usp.R[0];
        int thetaR = sVector[s].Usp.R[1];
        int aR = sVector[s].Uip.R[0];
        HANDLE_ERROR(cudaMallocArray(&arrayUsp, &channelDescUsp, dR, thetaR));
        HANDLE_ERROR(cudaMallocArray(&arrayUip, &channelDescUip, aR, thetaR));

        texUsp.addressMode[0] = cudaAddressModeMirror;
        texUsp.addressMode[1] = cudaAddressModeMirror;
        texUsp.filterMode     = cudaFilterModeLinear;
        texUsp.normalized     = false;

        texUip.addressMode[0] = cudaAddressModeMirror;
        texUip.addressMode[1] = cudaAddressModeMirror;
        texUip.filterMode     = cudaFilterModeLinear;
        texUip.normalized     = false;
        HANDLE_ERROR(cudaBindTextureToArray(texUsp, arrayUsp, channelDescUsp));
        HANDLE_ERROR(cudaBindTextureToArray(texUip, arrayUip, channelDescUip));

        //copy the LUT to the Usp texture
        HANDLE_ERROR( cudaMemcpy2DToArray(arrayUsp, 0, 0, sVector[s].Usp.x_hat, dR*sizeof(float2), dR*sizeof(float2), thetaR, cudaMemcpyDeviceToDevice));
        HANDLE_ERROR( cudaMemcpy2DToArray(arrayUip, 0, 0, sVector[s].Uip.x_hat, aR*sizeof(float2), aR*sizeof(float2), thetaR, cudaMemcpyDeviceToDevice));

        gpuScalarUpLut<<<dimGrid, dimBlock>>>(U.x_hat,
                                              gpuk,
                                              nWaves,
                                              2 * PI / lambda,
                                              sVector[s].a,
                                              sVector[s].d_min,
                                              sVector[s].d_max,
                                              focus,
                                              sVector[s].p,
                                              subA,
                                              pos,
                                              U.R[0],
                                              U.R[1],
                                              dR,
                                              aR,
                                              thetaR);

        //        bsVector c = sVector[s].p - focus;

        //        ptype *d_d, *d_thetai, *d_ai, *d_di;
        //        bsComplex *d_scaledPhase;
        //        HANDLE_ERROR(cudaMalloc( (void**)&d_thetai, sizeof(ptype) * U.R[0]*U.R[1]*inWaves.size()));

        //        HANDLE_ERROR(cudaMalloc( (void**)&d_d, sizeof(ptype) * U.R[0]*U.R[1]));
        //        HANDLE_ERROR(cudaMalloc( (void**)&d_ai, sizeof(ptype) * U.R[0]*U.R[1]));
        //        HANDLE_ERROR(cudaMalloc( (void**)&d_di, sizeof(ptype) * U.R[0]*U.R[1]));

        //        int sizeScaledPhase=(int) ceil(float(nWaves)/(SQRT_BLOCK*SQRT_BLOCK))*SQRT_BLOCK*SQRT_BLOCK;

        //        HANDLE_ERROR(cudaMalloc( (void**)&d_scaledPhase, sizeof(bsComplex) * sizeScaledPhase));

        //        ptype u = 1;
        //        ptype v = 1;
        //        //get the rtsPoint in world space and then the r vector
        //        bsPoint p = pos(u, v);
        //        bsVector r = p - sVector[s].p;

        //        std::cout<<"p: "<<p.toStr()<<std::endl;
        //        std::cout<<"ps: "<<sVector[s].p.toStr()<<std::endl;
        //        std::cout<<"r: "<<r.toStr()<<std::endl;


        //        precomputeLUT<<<dimGrid, dimBlock>>>(dR, aR,  sVector[s].d_max,  sVector[s].d_min, pos, nWaves,
        //                                             d_thetai, gpuk,
        //                                             d_ai, d_di, d_d, sVector[s].p, U.R[0], U.R[1], thetaR,
        //                sVector[s].a, d_scaledPhase, 2 * PI / lambda, c, subA);

        //        //   printf("\t\t shared memory size %d\n", sizeScaledPhase);

        //        //sizeof(bsComplex)*nWaves
        //        optPc_gpuScalarUpLut<<<dimGrid, dimBlock , sizeof(bsComplex)*sizeScaledPhase>>>(d_d,
        //                                                                                        U.x_hat,
        //                                                                                        nWaves,
        //                                                                                        sVector[s].a,
        //                                                                                        U.R[0],
        //                                                                                        U.R[1],
        //                                                                                        d_ai,
        //                                                                                        d_di,
        //                                                                                        d_thetai,
        //                                                                                        d_scaledPhase);

        cudaFreeArray(arrayUsp);
        cudaFreeArray(arrayUip);
        //cudaFree(d_thetai);
//        /cudaFree(d_d);
//        cudaFree(d_ai);
//        cudaFree(d_di);
//        cudaFree(d_scaledPhase);

    }


    //store the time to compute the scattered field
    t_Us = gpuStopTimer();

    //free monte-carlo samples
    cudaFree(gpuk);

}



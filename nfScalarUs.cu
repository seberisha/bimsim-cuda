#include "nearfield.h"
#include "rts/math/spherical_bessel.h"
#include "rts/math/legendre.h"
#include <stdlib.h>
#include "rts/cuda/error.h"
#include "rts/cuda/timer.h"

#define THREADS_PER_BLOCK          SQRT_BLOCK*SQRT_BLOCK
#if __CUDA_ARCH__ >= 200
    #define MY_KERNEL_MAX_THREADS  (2 * THREADS_PER_BLOCK)
    #define MY_KERNEL_MIN_BLOCKS   16
#else
    #define MY_KERNEL_MAX_THREADS  THREADS_PER_BLOCK
    #define MY_KERNEL_MIN_BLOCKS   16
#endif

__global__ void precompute(bsComplex n, bsRect ABCD, int nk, ptype *kd, bsComplex *knd, ptype *cos_theta, ptype kmag, bsVector* k, ptype *d, bsPoint ps, int uR, int vR)
{
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
    d[i] = r.len();
    knd[i] = kmag * d[i] * n;
    kd[i] = kmag*d[i];

    //normalize the direction vectors and find their inner product
    r = r.norm();

    for(int iw = 0; iw < nk; iw++)
    {
        cos_theta[iw*uR*vR + i] = k[iw].dot(r);
    }

}

__device__ bsComplex calc_Us(ptype kd, ptype cos_theta, int Nl, bsComplex* B)
{
	//initialize the spherical Bessel functions
	ptype j[2];
	rts::init_sbesselj<ptype>(kd, j);
	ptype y[2];
	rts::init_sbessely<ptype>(kd, y);

	//initialize the Legendre polynomial
	ptype P[2];
	rts::init_legendre<ptype>(cos_theta, P[0], P[1]);

	//initialize the spherical Hankel function
	bsComplex h((ptype)0, (ptype)0);

	//initialize the result
	bsComplex Us((ptype)0, (ptype)0);

	//for each order up to Nl
	for(int l=0; l<=Nl; l++)
	{
		if(l == 0)
		{
			h.r = j[0];
			h.i = y[0];
			Us += B[0] * h * P[0];
		}
		else
		{
			//shift the bessel functions and legendre polynomials
			if(l > 1)
			{
				rts::shift_legendre<ptype>(l, cos_theta, P[0], P[1]);
				rts::shift_sbesselj<ptype>(l, kd, j);
				rts::shift_sbessely<ptype>(l, kd, y);
			}

			h.r = j[1];
			h.i = y[1];
			Us += B[l] * h * P[1];


		}
	}
	return Us;
}


__device__ bsComplex opt_calc_Us(ptype kd, ptype cos_theta, int Nl, bsComplex* B)
{
    //initialize the spherical Bessel functions
    ptype j[2];
    //compute the first 2 bessel functions
    j[0] = sin(kd) / kd;

    j[1] = j[0] / kd - cos(kd) / kd;


    ptype y[2];

    //compute the first 2 bessel functions
    y[0] = -cos(kd) / kd;

    y[1] = y[0] / kd - sin(kd) / kd;

    //initialize the Legendre polynomial
    ptype P[2];
    //compute the first two Legendre polynomials
    P[0] = 1;
    P[1] = cos_theta;

    //initialize the spherical Hankel function
    bsComplex h((ptype)0, (ptype)0);

    //initialize the result
    bsComplex Us((ptype)0, (ptype)0);

    h.r = j[0];
    h.i = y[0];
    Us += B[0] * h * P[0];

    if(Nl==0)
        return;

    h.r = j[1];
    h.i = y[1];
    Us += B[1] * h * P[1];

    if(Nl==1)
        return;

    ptype temp=0;

    //for each order up to Nl
    for(int l=2; l<=Nl; l++)
    {
        //shift the bessel functions and legendre polynomials

       // rts::shift_legendre<ptype>(l, cos_theta, P[0], P[1]);

        //compute the next (order n) Legendre polynomial
        temp = ( (2 * l - 1) * cos_theta * P[1] - (l-1) * P[0] ) / l;

        //shift and add the new value to the array
        P[0] = P[1];
        P[1] = temp;


        //rts::shift_sbesselj<ptype>(l, kd, j);

        //compute the next (order n) Bessel function
        temp = ((2 * l - 1) * j[1])/kd - j[0];

        //if(n > stability*x)
        if(l > real(kd))
            if(real(temp) < RTS_BESSEL_CONVERGENCE_MIN || real(temp) > RTS_BESSEL_CONVERGENCE_MAX)
                temp = 0.0;

        //shift and add the new value to the array
        j[0] = j[1];
        j[1] = temp;


      //  rts::shift_sbessely<ptype>(l, kd, y);


        //compute the next (order n) Bessel function
        temp = ((2 * l - 1) * y[1])/kd - y[0];

        if(temp < RTS_BESSEL_MAXIMUM_FLOAT ||
           (l > kd && temp > 0))
        {
            temp = 0;
            y[1] = 0;
        }


        //shift and add the new value to the array
        y[0] = y[1];
        y[1] = temp;



        h.r = j[1];
        h.i = y[1];
        Us += B[l] * h * P[1];
    }
    return Us;
}


__device__ bsComplex calc_Ui(bsComplex knd, ptype cos_theta, int Nl, bsComplex* A)
{
	//calculate the internal field of a sphere

	bsComplex Ui((ptype)0, (ptype)0);

	//deal with rtsPoints near zero
	if(real(knd) < EPSILON_FLOAT)
	{
		//for(int l=0; l<Nl; l++)
		Ui = A[0];
        return Ui;
    }

	//initialize the spherical Bessel functions
	bsComplex j[2];
	rts::init_sbesselj<bsComplex>(knd, j);

	//initialize the Legendre polynomial
	ptype P[2];
	rts::init_legendre<ptype>(cos_theta, P[0], P[1]);

	//for each order up to Nl
	for(int l=0; l<=Nl; l++)
	{
		if(l == 0)
		{
			Ui += A[0] * j[0] * P[0];
		}
		else
		{
			//shift the bessel functions and legendre polynomials
			if(l > 1)
			{
				rts::shift_legendre<ptype>(l, cos_theta, P[0], P[1]);
				rts::shift_sbesselj<bsComplex>(l, knd, j);
			}

			Ui += A[l] * j[1] * P[1];


		}
	}
	return Ui;
}

__device__ bsComplex opt_calc_Ui(bsComplex knd, ptype cos_theta, int Nl, bsComplex* A)
{
    //calculate the internal field of a sphere

    bsComplex Ui((ptype)0, (ptype)0);

    //deal with rtsPoints near zero
    if(real(knd) < EPSILON_FLOAT)
    {
        //for(int l=0; l<Nl; l++)
        Ui = A[0];
        return Ui;
    }

    //initialize the spherical Bessel functions
    bsComplex j[2];
   // rts::init_sbesselj<bsComplex>(knd, j);

    //compute the first 2 bessel functions
    j[0] = sin(knd) / knd;

    j[1] = j[0] / knd - cos(knd) / knd;



    //initialize the Legendre polynomial
    ptype P[2];
    //rts::init_legendre<ptype>(cos_theta, P[0], P[1]);

    //compute the first two Legendre polynomials
    P[0] = 1;
    P[1] = cos_theta;


    Ui += A[0] * j[0] * P[0];

    if(Nl==0)
        return;

    Ui += A[1] * j[1] * P[1];

    if(Nl==1)
        return;

    ptype temp = 0;
    bsComplex ctemp(0,0);
    //for each order up to Nl
    for(int l=2; l<=Nl; l++)
    {

        //shift the bessel functions and legendre polynomials

        //rts::shift_legendre<ptype>(l, cos_theta, P[0], P[1]);

        //compute the next (order n) Legendre polynomial
        temp = ( (2 * l - 1) * cos_theta * P[1] - (l-1) * P[0] ) / l;

        //shift and add the new value to the array
        P[0] = P[1];
        P[1] = temp;


      //  rts::shift_sbesselj<bsComplex>(l, knd, j);



        //compute the next (order n) Bessel function
        ctemp = ((2 * l - 1) * j[1])/knd - j[0];

        //if(n > stability*x)
        if(l > real(knd))
            if(real(ctemp) < RTS_BESSEL_CONVERGENCE_MIN || real(ctemp) > RTS_BESSEL_CONVERGENCE_MAX)
                ctemp = 0.0;

        //shift and add the new value to the array
        j[0] = j[1];
        j[1] = ctemp;


        Ui += A[l] * j[1] * P[1];

    }
    return Ui;
}


__global__ void
//__launch_bounds__(MY_KERNEL_MAX_THREADS,MY_KERNEL_MIN_BLOCKS)
gpuScalarUsp(bsComplex* Us, bsVector* k, int nk, ptype kmag, bsPoint f, bsPoint ps, ptype a, bsComplex n, bsComplex* Beta, bsComplex* Alpha, int Nl, ptype A, bsRect ABCD, int uR, int vR)
{

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

    bsComplex sumUs(0, 0);
    //for each plane wave in the wave list

    for(int iw = 0; iw < nk; iw++)
    {
        //normalize the direction vectors and find their inner product
        r = r.norm();
        ptype cos_theta = k[iw].dot(r);

        //compute the phase factor for spheres that are not at the origin
        bsVector c = ps - f;


        bsComplex phase = exp(bsComplex(0, kmag * k[iw].dot(c)));
        //cout<<"\t\t green "<<green.abs()<<endl;
        //compute the internal field if we are inside a sphere
        if(d <= a)
        {
            bsComplex knd = kmag * d * n;
            sumUs += (1.0f/nk) * A * phase * calc_Ui(knd, cos_theta, Nl, Alpha);
        }
        //otherwise compute the scattered field
        else
        {
            //compute the argument for the spherical Hankel function
            ptype kd = kmag * d;
            sumUs += (1.0f/nk) * A * phase * calc_Us(kd, cos_theta, Nl, Beta);
        }

    }

    Us[i] += sumUs;


}


__global__ void
//__launch_bounds__(MY_KERNEL_MAX_THREADS,MY_KERNEL_MIN_BLOCKS)
opt_gpuScalarUsp(bsComplex* Us, bsVector* k, int nk, ptype kmag, bsVector c, bsPoint ps, ptype a, bsComplex n, bsComplex* Beta, bsComplex* Alpha, int Nl, ptype A, bsRect ABCD, int uR, int vR)
{

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

    bsComplex sumUs(0, 0);
    //for each plane wave in the wave list

    //normalize the direction vectors and find their inner product
    r = r.norm();
    //compute the internal field if we are inside a sphere
    if(d <= a)
    {
        for(int iw = 0; iw < nk; iw++)
        {
            sumUs += (1.0f/nk) * A * exp(bsComplex(0, kmag * k[iw].dot(c))) * opt_calc_Ui(kmag * d * n, k[iw].dot(r), Nl, Alpha);
        }
    }
    //otherwise compute the scattered field
    else
    {
        for(int iw = 0; iw < nk; iw++)
        {
            sumUs += (1.0f/nk) * A * exp(bsComplex(0, kmag * k[iw].dot(c))) * opt_calc_Us(kmag*d, k[iw].dot(r), Nl, Beta);

        }
    }
    Us[i] += sumUs;


}



__global__ void
//__launch_bounds__(MY_KERNEL_MAX_THREADS,MY_KERNEL_MIN_BLOCKS)
optPc_gpuScalarUsp(bsComplex* Us, bsVector* k, int nk, ptype kmag, bsVector c, bsComplex *knd, ptype* kd, ptype *d, ptype a, ptype *cos_theta, bsComplex* Beta, bsComplex* Alpha, int Nl, ptype A, int uR, int vR)
{

    //get the current coordinate in the plane slice
    int iu = blockIdx.x * blockDim.x + threadIdx.x;
    int iv = blockIdx.y * blockDim.y + threadIdx.y;

    //make sure that the thread indices are in-bounds
    if(iu >= uR || iv >= vR) return;

    //compute the index (easier access to the scalar field array)
    int i = iv*uR + iu;

    bsComplex sumUs(0, 0);
    //for each plane wave in the wave list
    ptype scale = (1.0f/nk);
    //compute the internal field if we are inside a sphere
    if(d[i] <= a)
    {
        for(int iw = 0; iw < nk; iw++)
        {
          //  bsComplex phase = exp(bsComplex(0, kmag * k[iw].dot(c)));
          //  sumUs+=0;
            sumUs += scale * A * exp(bsComplex(0, kmag * k[iw].dot(c))) * calc_Ui(knd[i], cos_theta[iw*uR*vR+i], Nl, Alpha);
                    //opt_calc_Ui(knd[i], cos_theta[iw*uR*vR+i], Nl, Alpha);
        }
    }
    //otherwise compute the scattered field
    else
    {
        for(int iw = 0; iw < nk; iw++)
        {
           // bsComplex phase = exp(bsComplex(0, kmag * k[iw].dot(c)));

            sumUs += scale * A * exp(bsComplex(0, kmag * k[iw].dot(c))) * calc_Us(kd[i],cos_theta[iw*uR*vR+i], Nl, Beta);
                    //opt_calc_Us(kd[i],cos_theta[iw*uR*vR+i], Nl, Beta);

        }
    }
    Us[i] += sumUs;


}


void nearfieldStruct::scalarUs()
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

//    cudaGetDeviceProperties(&deviceProp, device);
//    int threadsPerBlock =
//              (deviceProp.major >= 2 ?
//                        2 * THREADS_PER_BLOCK : THREADS_PER_BLOCK);

	//create one thread for each pixel of the field slice
	dim3 dimBlock(SQRT_BLOCK, SQRT_BLOCK);
    //printf("\t\t U.R: %d,%d\n",U.R[0],U.R[1]);
	dim3 dimGrid((U.R[0] + SQRT_BLOCK -1)/SQRT_BLOCK, (U.R[1] + SQRT_BLOCK - 1)/SQRT_BLOCK);

	//for each sphere
	int Nl;
	for(int s = 0; s<nSpheres; s++)
	{
		//get the number of orders for this sphere
		Nl = sVector[s].Nl;

		//allocate memory for the scattering coefficients
		bsComplex* gpuB;
		HANDLE_ERROR(cudaMalloc((void**) &gpuB, (Nl+1) * sizeof(bsComplex)));

		bsComplex* gpuA;
		HANDLE_ERROR(cudaMalloc((void**) &gpuA, (Nl+1) * sizeof(bsComplex)));

		//copy the scattering coefficients to the GPU
		HANDLE_ERROR(cudaMemcpy(gpuB, &sVector[s].B[0], (Nl+1) * sizeof(bsComplex), cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(gpuA, &sVector[s].A[0], (Nl+1) * sizeof(bsComplex), cudaMemcpyHostToDevice));

		//if we are computing a plane wave, call the gpuScalarUfp function
		sphere S = sVector[s];
		bsVector* gpuk;

		if(planeWave)
		{
           // std::cout<<"void nearfieldStruct::scalarUs(): single plane wave"<<std::endl;
            //if this is a single plane wave, assume it goes along direction k (copy the k vector to the GPU)
            HANDLE_ERROR(cudaMalloc( (void**)&gpuk, sizeof(bsVector) ));
            HANDLE_ERROR(cudaMemcpy( gpuk, &k, sizeof(bsVector), cudaMemcpyHostToDevice));
			gpuScalarUsp<<<dimGrid, dimBlock>>>(U.x_hat,
												gpuk,
												1,
												2 * PI / lambda,
												focus,
												sVector[s].p,
												sVector[s].a,
												sVector[s].n,
												gpuB,
												gpuA,
												Nl,
												A,
												pos,
												U.R[0],
												U.R[1]);
            HANDLE_ERROR(cudaFree(gpuk));
		}
		//otherwise copy all of the monte-carlo samples to the GPU and compute
		else
		{
//            std::cout<<"void nearfieldStruct::scalarUs(): multiple plane wave"<<std::endl;
            HANDLE_ERROR(cudaMalloc( (void**)&gpuk, sizeof(bsVector) * inWaves.size() ));
            HANDLE_ERROR(cudaMemcpy( gpuk, &inWaves[0], sizeof(bsVector) * inWaves.size(), cudaMemcpyHostToDevice));

            //compute the amplitude that makes it through the condenser
            ptype subA = 2 * PI * A * ( (1 - cos(asin(condenser[1]))) - (1 - cos(asin(condenser[0]))) );

//            gpuScalarUsp<<<dimGrid, dimBlock>>>(U.x_hat,
//												gpuk,
//												inWaves.size(),
//												2 * PI / lambda,
//												focus,
//												sVector[s].p,
//												sVector[s].a,
//												sVector[s].n,
//												gpuB,
//												gpuA,
//												Nl,
//												subA,
//												pos,
//												U.R[0],
//												U.R[1]);
            bsVector c = sVector[s].p - focus;

            ptype *d_kd, *d_cos_theta, *d_d;
            bsComplex *d_knd;
            HANDLE_ERROR(cudaMalloc( (void**)&d_kd, sizeof(ptype) * U.R[0]*U.R[1]));
            HANDLE_ERROR(cudaMalloc( (void**)&d_knd, sizeof(bsComplex) * U.R[0]*U.R[1]));
            HANDLE_ERROR(cudaMalloc( (void**)&d_cos_theta, sizeof(ptype) * U.R[0]*U.R[1]*inWaves.size()));

            HANDLE_ERROR(cudaMalloc( (void**)&d_d, sizeof(ptype) * U.R[0]*U.R[1]));




            precompute<<<dimGrid, dimBlock>>>(sVector[s].n,
                                              pos,
                                              inWaves.size(),
                                              d_kd,
                                              d_knd,
                                              d_cos_theta,
                                              2 * PI / lambda,
                                              gpuk,
                                              d_d,
                                              sVector[s].p,
                                              U.R[0],
                                              U.R[1]);




            //HANDLE_ERROR(cudaMemcpy( gpuk, &inWaves[0], sizeof(bsVector) * inWaves.size(), cudaMemcpyHostToDevice));



            optPc_gpuScalarUsp<<<dimGrid, dimBlock>>>(U.x_hat,
                                                gpuk,
                                                inWaves.size(),
                                                2 * PI / lambda,
                                                c,
                                                d_knd,
                                                d_kd,
                                                d_d,
                                                sVector[s].a,
                                                d_cos_theta,
                                                gpuB,
                                                gpuA,
                                                Nl,
                                                subA,
                                                U.R[0],
                                                U.R[1]);

            HANDLE_ERROR(cudaFree(gpuk));
            HANDLE_ERROR(cudaFree(d_kd));
            HANDLE_ERROR(cudaFree(d_knd));
            HANDLE_ERROR(cudaFree(d_cos_theta));
            HANDLE_ERROR(cudaFree(d_d));


		}

		//free memory for scattering coefficients
        HANDLE_ERROR(cudaFree(gpuA));
        HANDLE_ERROR(cudaFree(gpuB));
	}


    //store the time to compute the scattered field
	t_Us = gpuStopTimer();


}

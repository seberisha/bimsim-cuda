#include "fieldslice.h"
#include "dataTypes.h"

#include "cufft.h"

#include <iostream>
#include <stdlib.h>
using namespace std;

fieldslice::fieldslice(unsigned int x_size, unsigned int y_size)
{
    x_hat = y_hat = z_hat = NULL;

	//save the slice resolution
	R[0] = x_size;
	R[1] = x_size;

	scalarField = true;

	init_gpu();


}

void fieldslice::toAngularSpectrum()
{
    cufftHandle plan;
	cufftResult result;

    //create a CUFFT plan handle
#ifdef PRECISION_SINGLE
    result = cufftPlan2d(&plan, R[0], R[1], CUFFT_C2C);
#elif defined PRECISION_DOUBLE
	 result = cufftPlan2d(&plan, R[0], R[1], CUFFT_Z2Z);
#endif

	if(result != CUFFT_SUCCESS)
    {
        cout<<"Error creating CUFFT plan for computing the angular spectrum."<<endl;
        exit(1);
    }

#ifdef PRECISION_SINGLE
    result = cufftExecC2C(plan, (cufftComplex*)x_hat, (cufftComplex*)x_hat, CUFFT_FORWARD);
#elif defined PRECISION_DOUBLE
    result = cufftExecZ2Z(plan, (cufftDoubleComplex*)x_hat, (cufftDoubleComplex*)x_hat, CUFFT_FORWARD);
#endif
	if(result != CUFFT_SUCCESS)
    {
        cout<<"Error executing the CUFFT forward FFT to compute the angular spectrum."<<endl;
        exit(1);

    }

    cufftDestroy(plan);

}

void fieldslice::fromAngularSpectrum()
{
    cufftHandle plan;
	cufftResult result;

    //create a CUFFT plan handle
#ifdef PRECISION_SINGLE
    result = cufftPlan2d(&plan, R[0], R[1], CUFFT_C2C);
#elif defined PRECISION_DOUBLE
	 result = cufftPlan2d(&plan, R[0], R[1], CUFFT_Z2Z);
#endif

	if(result != CUFFT_SUCCESS)
    {
        cout<<"Error creating CUFFT plan for computing the angular spectrum."<<endl;
        exit(1);
    }

#ifdef PRECISION_SINGLE
    result = cufftExecC2C(plan, (cufftComplex*)x_hat, (cufftComplex*)x_hat, CUFFT_INVERSE);
#elif defined PRECISION_DOUBLE
    result = cufftExecZ2Z(plan, (cufftDoubleComplex*)x_hat, (cufftDoubleComplex*)x_hat, CUFFT_INVERSE);
#endif
	if(result != CUFFT_SUCCESS)
    {
        cout<<"Error executing the CUFFT forward FFT to compute the angular spectrum."<<endl;
        exit(1);

    }

    //divide the field by the number of values in the field
    ScaleField( 1.0 / (R[0] * R[1]));

    cufftDestroy(plan);

}

fieldslice::fieldslice()
{
	R[0] = R[1] = 0;
	x_hat = y_hat = z_hat = NULL;
}



fieldslice::~fieldslice()
{
	kill_gpu();
}

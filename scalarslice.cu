#include "scalarslice.h"

#include "rts/cuda/error.h"
#include "cublas_v2.h"
#include "rts/envi/envi.h"

scalarslice::scalarslice(int x, int y)
{
	//set the resolution
	R[0] = x;
	R[1] = y;

	//allocate memory on the GPU
	HANDLE_ERROR(cudaMalloc( (void**)&S, sizeof(ptype) * x * y ));

    //std::cout<<"Scalerslice created."<<std::endl;
}



scalarslice::scalarslice(int x, int y, int psSamples)
{
    //set the resolution
    R[0] = x;
    R[1] = y;

    //allocate memory on the GPU
    HANDLE_ERROR(cudaMalloc( (void**)&S, sizeof(ptype) * x * y ));

    //allocate memoru on the GPU for the radius of the detector image
    HANDLE_ERROR(cudaMalloc((void**)&line, sizeof(ptype)*psSamples));

    //std::cout<<"Scalerslice created."<<std::endl;
}


scalarslice::scalarslice()
{
    R[0] = R[1] = 0;
    S = NULL;

    //std::cout<<"Scalerslice created (default)."<<std::endl;
}

scalarslice::scalarslice(ptype s)
{
    R[0] = R[1] = 1;
    *S = s;
    //std::cout<<"Scalerslice created (default)."<<std::endl;
}




scalarslice::~scalarslice()
{
	if(S != NULL)
		HANDLE_ERROR(cudaFree(S));
	S = NULL;
	R[0] = R[1] = 0;

	//std::cout<<"Scalerslice destroyed."<<std::endl;
}

//assignment operator
scalarslice & scalarslice::operator= (const scalarslice & rhs)
{
    //de-allocate any existing GPU memory
    if(S != NULL)
		HANDLE_ERROR(cudaFree(S));

    //copy the slice resolution
    R[0] = rhs.R[0];
    R[1] = rhs.R[1];

    //allocate the necessary memory
    HANDLE_ERROR(cudaMalloc(&S, sizeof(ptype) * R[0] * R[1]));

    //copy the slice
    HANDLE_ERROR(cudaMemcpy(S, rhs.S, sizeof(ptype) * R[0] * R[1], cudaMemcpyDeviceToDevice));


    std::cout<<"Assignment operator."<<std::endl;

	return *this;

}

void scalarslice::toImage(std::string filename, ptype vmin, ptype vmax, stim::colormapType cmap)
{
	stim::gpu2image<ptype>(S, filename, R[0], R[1], vmin, vmax, cmap);
}

void scalarslice::toImage(std::string filename, bool positive, stim::colormapType cmap)
{
   // std::cout<<"toImage(std::string filename, bool positive, rts::colormapType cmap)"<<std::endl;

    cublasStatus_t stat;
    cublasHandle_t handle;

    //create a CUBLAS handle
    stat = cublasCreate(&handle);
    if(stat != CUBLAS_STATUS_SUCCESS)
    {
        std::cout<<"CUBLAS Error: initialization failed"<<std::endl;
        exit(1);
    }

    //find the index of the value with maximum magnitude
    int N = R[0] * R[1];
    int result;
#ifdef PRECISION_SINGLE
    stat = cublasIsamax(handle, N, S, 1, &result);
#elif defined PRECISION_DOUBLE
	stat = cublasIdamax(handle, N, S, 1, &result);
#endif


	//adjust for 1-based indexing
	result -= 1;

	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		std::cout<<"CUBLAS Error: failure finding maximum value."<<std::endl;
		exit(1);
	}



    //retrieve the maximum value
    ptype maxVal;
    HANDLE_ERROR(cudaMemcpy(&maxVal, S + result, sizeof(ptype), cudaMemcpyDeviceToHost));

    //destroy the CUBLAS handle
    cublasDestroy(handle);

    std::cout<<"\t\t positive, result, maxVal: "<<positive<<" , "<<result<<" , "<<" , "<<maxVal<<endl;

    //output the image
    if(positive)
        toImage(filename, 0, maxVal, cmap);
    else
        toImage(filename, -abs(maxVal), abs(maxVal), cmap);
}

void scalarslice::toEnvi(std::string filename, ptype wavelength, bool append)
{
    std::string mode;
    if(append) mode = "a";
    else       mode = "w";

    //open the ENVI file
    EnviFile outfile(filename, mode);

    //get the scalar slice from the GPU to the CPU
    int memsize = sizeof(ptype) * R[0] * R[1];
    ptype* cpuData = (ptype*) malloc( memsize );
    HANDLE_ERROR(cudaMemcpy( cpuData, S, memsize, cudaMemcpyDeviceToHost));

    //add a band to the ENVI file
    outfile.addBand(cpuData, R[0], R[1], wavelength);

    outfile.close();


}

void scalarslice::clear()
{
	//this function sets the slice to zero
	if(S != NULL)
    {
		//HANDLE_ERROR(cudaMalloc( (void**)&S, sizeof(ptype) * R[0] * R[1] ));
		HANDLE_ERROR(cudaMemset(S, 0, sizeof(ptype) * R[0] * R[1]));
        //cout<<"\t\t S memset"<<endl;
    }

//    for (int i=0; i<7;++i)
//        cout<<S[i]<<endl;
}

void scalarslice::clearLine(int radiusLength)
{
    //this function sets the radius of the slice to zero
    if(line != NULL)
    {
        //HANDLE_ERROR(cudaMalloc( (void**)&S, sizeof(ptype) * R[0] * R[1] ));
        HANDLE_ERROR(cudaMemset(line, 0, sizeof(ptype) * radiusLength));
    }
}


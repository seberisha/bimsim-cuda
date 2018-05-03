#include "nearfield.h"
#include <time.h>
#include <math.h>

#ifdef _WIN32
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#endif

int bessjyv_sph(int v, double z, double &vm, double* cjv,
    double* cyv, double* cjvp, double* cyvp);

nearfieldStruct::nearfieldStruct()
{
    scalarSim = true;
	planeWave = false;

    //sb: changing this to false--- lut_us = true;
    lut_us = true;
	lut_uf = false;

	nWaves = 0;
}

void nearfieldStruct::init()
{
	//set the field parameters
	U.scalarField = scalarSim;
	Uf.scalarField = scalarSim;

	//initialize dynamic memory
	U.init_gpu();
	Uf.init_gpu();
}

void nearfieldStruct::destroy()
{
	U.kill_gpu();
	Uf.kill_gpu();
}

void nearfieldStruct::setPos(bsPoint pMin, bsPoint pMax, bsVector normal)
{
    std::cout<<"\t\t pMin, pMax, normal "<<pMin.toStr()<<" , "<<pMax.toStr()<<" , "<<normal.toStr()<<std::endl;
	pos = rts::quad<ptype, 3>(pMin, pMax, normal);
}

void nearfieldStruct::setRes(int x_res, int y_res)
{
	U.R[0] = Uf.R[0] = x_res;
	U.R[1] = Uf.R[1] = y_res;
}

std::string nearfieldStruct::toStr()
{
	std::stringstream ss;

    ss<<"---------Timings-------------"<<std::endl;
    ss<<"Uf = "<<t_Uf<<"ms"<<std::endl;
    ss<<"Us = "<<t_Us<<"ms"<<std::endl;
    ss<<std::endl;

    ss<<"------Field Parameters-------"<<std::endl;
    ss<<"Wavelength: "<<lambda<<"um"<<std::endl;
    ss<<"K Vector (x, y, z): "<<k<<std::endl;
    ss<<"K Vector (r, theta, phi): "<<k.cart2sph()<<std::endl;
    ss<<"Condenser NA: "<<condenser[0]<<" to "<<condenser[1]<<std::endl;
    ss<<"Number of points and rings: "<<extSource[0]<<" to "<<extSource[1]<<std::endl;
    ss<<"Focal Point: "<<focus[0]<<", "<<focus[1]<<", "<<focus[2]<<std::endl;
    //ss<<"Field Slice: "<<std::endl;
    if(lut_us)
        ss<<"LUT Parameters --- min: "<<d_min<<"   max: "<<d_max<<std::endl;
    ss<<pos<<std::endl;

    ss<<std::endl<<"---------Materials-----------"<<std::endl;
    ss<<"Number of Materials: "<<mVector.size()<<std::endl;
    ss<<"Refractive Indices at lambda = "<<lambda<<"um"<<std::endl;
    //output each material
    for(unsigned int m=0; m<mVector.size(); m++)
        ss<<" "<<m<<": "<<mVector[m](lambda)<<std::endl;

    ss<<"---------Spheres-------------"<<std::endl;
    ss<<"Number of Spheres: "<<sVector.size()<<std::endl;
    //output each sphere
    for(unsigned int s=0; s<sVector.size(); s++)
        ss<<sVector[s].toStr()<<std::endl;


    //output only lambda, n, k
    ss<<lambda<<","<<mVector[0](lambda)<<std::endl;

	return ss.str();
}

//generate monte-carlo waves
void nearfieldStruct::calcWaves()
{
    inWaves.resize(nWaves);

    //re-seed the random number generator
    //srand(time(NULL));
	srand(NULL);

    //calculate the monte-carlo samples
    mcSampleNA(&inWaves[0], nWaves, k, condenser[0], condenser[1]);


}


void nearfieldStruct::calcSpheres()
{
    //calculate all of the constants necessary to evaluate the scattered field
	//estimate the order required to represent the scattered field for each sphere



	//for each sphere
	for(unsigned int i=0; i<sVector.size(); i++)
	{
		//a = sVector[i].a;

		//calculate the required order
		sVector[i].calcNl(lambda);


		//set the refractive index for the sphere
		int imat = sVector[i].iMaterial;
        rts::complex<ptype> n = mVector[imat](lambda);

		//calculate the scattering coefficients
		sVector[i].calcCoeff(lambda, n);

		//save the refractive index
		sVector[i].n = n;

		//if the LUT is used, calculate Usp(theta, r)
		if(lut_us)
		{
			sVector[i].calcUp(lambda, n, pos, std::max(U.R[0], U.R[1]));
        }


	}

}

void nearfieldStruct::calcUs()
{

    if(scalarSim)
    {
        if(lut_us)
            scalarUpLut();
        else
            scalarUs();

//        if(lut_us)
//            std::cout<<"Using LUT for Us"<<std::endl;
      //  std::cout<<"Calculate Us Scalar Sim."<<std::endl;
    }
    else
    {
        std::cout<<"Calculate Us Vector Sim."<<std::endl;
    }
}

void nearfieldStruct::calcUf()
{
    if(scalarSim)
    {
        if(lut_uf)
            scalarUfLut();
        else
            scalarUf();

        if(lut_uf)
            std::cout<<"Using LUT for Uf"<<std::endl;
       // std::cout<<"void nearfieldStruct::calcUf(): Calculate Uf Scalar Sim."<<std::endl;
    }
    else
    {
        std::cout<<"Calculating Uf Vector Sim..."<<std::endl;

        if(lut_uf)
            vectorUfLut();
        else
            vectorUf();


    }
}

void nearfieldStruct::Simulate()
{
    //initialize timings
    t_Uf = 0;
    t_Us = 0;
    //printf("\n\t\t\t Before calcWaves\n");
	//compute a set of plane waves for Monte-Carlo simulation
	calcWaves();
    //printf("\n\t\t\t Done calcWaves\n");

    //the near field has to be simulated no matter what the output rtsPoint is
    calcUf();
    calcSpheres();
    calcUs();
     
    //U.Mag().toImage("Us-nf.bmp");
    sumUf();
    /*
    if(verbose)
    {
        U.Mag().toImage("U-nf.bmp");
        Uf.Mag().toImage("Uf-nf.bmp");
     }
    */	
}



void nearfieldStruct::calcBesselLut(ptype* j, ptype d_min, ptype d_max, int dR)
{
    /*Compute the look-up-table for spherical bessel functions used for the incident field
        j    =   (Nl + 1) x aR array of values
        aR      =   resolution of j
    */

	//compute the wavenumber
	ptype k = 2 * PI / lambda;
	unsigned int Nl = m;

    //allocate space for the Bessel functions of the first and second kind (and derivatives -- which will be ignored)
    int bytes = sizeof(double) * (Nl + 1);
    double* cjv_kd = (double*)malloc(bytes);
    double* cyv_kd = (double*)malloc(bytes);
    double* cjvp_kd = (double*)malloc(bytes);
    double* cyvp_kd = (double*)malloc(bytes);

    //compute the bessel functions using the CPU-based algorithm
    double vm;

    //for each sample along r
    ptype dr = (d_max - d_min) / (dR - 1);
    ptype d;
    ptype jv;
    for(int id = 0; id < dR; id++)
    {
        d = id * dr + d_min;
        double kd = k*d;
        bessjyv_sph(Nl, kd, vm, cjv_kd, cyv_kd, cjvp_kd, cyvp_kd);

        //copy the double data to the bsComplex array
        for(unsigned int l=0; l<=Nl; l++)
		{
            jv = cjv_kd[l];
			if(std::isnan(jv) || std::isinf(jv))
			{
                if(kd == 0 && l == 0)
                    jv = 1;
                else
                    jv = 0;
            }
            j[id * (Nl+1) + l] = jv;
		}
    }

    if(verbose)
    {
        std::ofstream outfile("uf_besselout.txt");
        for(int ir = 0; ir < dR; ir++)
        {
            outfile<<ir*dr + d_min<<std::endl;
            for(unsigned int l = 0; l<=Nl; l++)
            {
                outfile<<j[ir * (Nl+1) + l]<<" -- ";
            }
            outfile<<std::endl;
        }
        outfile.close();
    }

}

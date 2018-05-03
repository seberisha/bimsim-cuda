#include "sphere.h"
#include "rts/optics/material.h"
#include <vector>

#include "nearfield.h"

#include "microscope.h"
microscopeStruct* SCOPE;

#include "fieldslice.h"

#include "fileout.h"
#include "arguments.h"
#include "rts/tools/arguments.h"
#include "montecarlo.h"
#include "rts/math/point.h"
#include "rts/math/spherical_bessel.h"
#include "rts/math/matrix.h"
#include "rts/math/quaternion.h"

#include "rts/envi/envi.h"

#include "warnings.h"

#include "planewave.h"

bool gui = false;

#ifdef _WIN32
bool ansi = false;
#else
bool ansi = true;
#endif

fileoutStruct gFileOut;
bool verbose = false;
using namespace std;

int cbessjyva(double v,complex<double> z,double &vm,complex<double>*cjv,
    complex<double>*cyv,complex<double>*cjvp,complex<double>*cyvp);

int main(int argc, char *argv[])
{


    cudaSetDevice(0);
    //arguments test
    rts::arglist args;
    SetArguments(args);

    //parse the input arguments 
    args.parse(argc, argv);

    SCOPE = new microscopeStruct(); 
    //load the user specified parameters into the simulation
    LoadParameters(args);


    //activate ansi output if specified
    args.set_ansi(ansi);

    //display help and exit
    if(args("help"))
    {
    	cout<<args.toStr()<<endl;
	exit(1);
    }

    //initialize GPU memory for fields
    SCOPE->init();

    gFileOut.Save(SCOPE);

    if(verbose)
    	OutputOptions();

    SCOPE->destroy();
}

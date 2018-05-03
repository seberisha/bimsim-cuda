#include "rts/tools/arguments.h"

#include "rts/optics/material.h"

#include "nearfield.h"
#include "microscope.h"
#include "stim/visualization/colormap.h"
#include "fileout.h"
extern microscopeStruct* SCOPE;
extern fileoutStruct gFileOut;

//default values
#include "defaults.h"

#include <string>
#include <sstream>
#include <fstream>
#include <limits>

extern bool verbose;
extern bool gui;

#ifdef _WIN32
	extern bool ansi;
#endif

void SetArguments(rts::arglist &args)
{
    args.section("Interface Flags");
	args.add("help", "prints this help");
	//args.add("gui", "run using the Qt GUI");
	args.add("verbose", "verbose output");

	args.section("Output Parameters");
	args.add("vector", "run a vector field simulation");
	args.add("intensity", "output measured intensity (filename)", DEFAULT_INTENSITY_FILE);
    args.add("incident", "output measured incident field (filename)", DEFAULT_INCIDENT_FILE);

	args.add("absorbance", "output measured absorbance (filename)", DEFAULT_ABSORBANCE_FILE);
    args.add("absorbance-spec", "output measured absorbance spectrum (filename)", DEFAULT_ABSORBANCE_SPECTRUM_FILE);

	args.add("transmittance", "output measured transmittance (filename)", DEFAULT_TRANSMITTANCE_FILE);
	args.add("far-field", "output far-field at detector (filename)", DEFAULT_FAR_FILE);
	args.add("near-field", "output field at focal plane (filename)", DEFAULT_NEAR_FILE);
	args.add("extended-source", "image of source at focus (filename)", DEFAULT_EXTENDED_SOURCE);

    args.add("interpolate", "option to interpolate the point sources", DEFAULT_INTERPOLATE  );


	args.add("output-type", "output field value", DEFAULT_FIELD_TYPE, "magnitude, polarization, real, imaginary");
	args.add("colormap", "colormap", DEFAULT_COLORMAP, "gray, brewer");
	args.add("append", "append result to an existing (binary) file");

	args.section("Sphere Parameters");
	args.add("spheres", "sphere position", "", "--spheres x y z a m");
    args.add("sphere-file", "sphere file:", "", "[x y z radius material]");
    args.add("materials", "refractive indices as n, k pairs", DEFAULT_MATERIAL, "--materials n0 k0 n1 k1 n2 k2");
    args.add("material-file", "material file", "", "[lambda n k]");

    args.section("Optics");
    args.add("lambda", "incident wavelength (micrometers)", DEFAULT_LAMBDA);
    args.add("nu", "incident frequency (in cm^-1)\n(if specified, lambda is ignored)");
    args.add("k", "k-vector direction in spherical coordinates", "", "--k theta phi; theta = [0 2*pi], phi = [0 pi]");
    args.add("amplitude", "incident field amplitude", DEFAULT_AMPLITUDE);
    args.add("condenser", "condenser numerical aperature\nA pair of values specify an inner obscuration", DEFAULT_CONDENSER);
    args.add("objective", "objective numerical aperature\nA pair of values specify an inner obscuration", DEFAULT_OBJECTIVE);
    args.add("focus", "focal position for the incident point source", DEFAULT_FOCUS);
    args.add("plane-wave", "simulates an incident plane wave");
    args.add("points-rings", "pair of values specying the number of rings and the number of point sources per ring", DEFAULT_POINTS_RINGS);

    args.section("Imaging Parameters");
    args.add("resolution", "resolution of the detector", DEFAULT_SLICE_RES);
	args.add("plane-lower-left", "lower-left position of the image plane", DEFAULT_PLANE_MIN);
	args.add("plane-upper-right", "upper-right position of the image plane", DEFAULT_PLANE_MAX);
	args.add("plane-normal", "normal for the image plane", DEFAULT_PLANE_NORM);
    args.add("xy", "specify an x-y axis-aligned image (standard microscope)");
    args.add("xz", "specify a x-z axis-aligned image");
    args.add("yz", "specify a y-z axis-aligned image");

    args.section("Sampling Parameters");
    args.add("samples", "Monte-Carlo samples used to compute Us", DEFAULT_SAMPLES);
	args.add("padding", "FFT padding for the objective bandpass", DEFAULT_PADDING);
	args.add("supersample", "super-sampling rate for the detector field", DEFAULT_SUPERSAMPLE);
	args.add("field-order", "order of the incident field", DEFAULT_FIELD_ORDER);
	args.add("seed", "seed for the Monte-Carlo random number generator");
	args.add("recursive", "evaluate all Bessel functions recursively");
	args.add("recursive-us", "evaluate scattered-field Bessel functions recursively");
	args.add("lut-uf", "evaluate the focused-field using a look-up table");

}

void lFlags(rts::arglist args)
{

    //flag for verbose output
	if(args("verbose"))
        verbose = true;

    if(args("recursive"))
    {
        SCOPE->nf.lut_us = false;
        SCOPE->nf.lut_uf = false;
    }
    else if(args("recursive-us"))
    {
        SCOPE->nf.lut_us = false;
    }
    else if(args("lut-uf"))
    {
        SCOPE->nf.lut_uf = true;
    }

}

void lWavelength(rts::arglist args)
{
    //load the wavelength
	if(args("nu"))
	{
		//wavelength is given in wavenumber - transform and flag
		SCOPE->nf.lambda = 10000/args["nu"].as_float();
		gFileOut.wavenumber = true;
	}
	//otherwise we are using lambda = wavelength
	else
	{
		SCOPE->nf.lambda = args["lambda"].as_float();
     //   std::cout<<"\t\t nf.lambda "<<SCOPE->nf.lambda<<std::endl;

		gFileOut.wavenumber = false;
	}
}

static void lMaterials(rts::arglist args)
{
    //SB:both if statements evaluate to true if given an argument of material-file

	//if file names are specified, load the materials
	if(args("material-file"))
	{
      //  std::cout<<"\t\t it's a material file"<<std::endl;
		rts::argument matfiles = args["material-file"];
      //  std::cout<<"\t\t matfiles "<<matfiles.as_text()<<std::endl;
		int nMats = matfiles.nargs();
     //   std::cout<<"\t\t nMats: "<<nMats<<std::endl;

        for(unsigned int i=0; i<nMats; i++)
        {
            //load the file into a string
            //std::ifstream ifs(filenames[i].c_str());

            //std::string instr((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

            //load the list of spheres from a string
            rts::material<ptype> newM(matfiles.as_text(i));

            SCOPE->nf.mVector.push_back(newM);
        }

    }else //if materials are specified at the command line if(args("materials"))
    {
       // std::cout<<"\t\t it's materials"<<std::endl;

        rts::argument mats = args["materials"];
        int nMats = mats.nargs();
      //  std::cout<<"\t\t nMats: "<<nMats<<std::endl;

        if(nMats == 1)
        {
            rts::material<ptype> newM(SCOPE->nf.lambda, mats.as_float(), 0);
            SCOPE->nf.mVector.push_back(newM);
        }
        else if(nMats %2 != 0)
        {
            cout<<"BIMSim Error: materials must be specified in n, k pairs"<<endl;
            exit(1);
        }
        else
        {
            for(unsigned int i=0; i<nMats; i+=2)
            {
                rts::material<ptype> newM(SCOPE->nf.lambda, mats.as_float(i), mats.as_float(i+1));
             //   std::cout<<"\t\t===>newM(SCOPE->nf.lambda, mats.as_float(i), mats.as_float(i+1) "<<
               //         SCOPE->nf.lambda <<" "<< mats.as_float(i) <<" " << mats.as_float(i+1)<<std::endl;
                SCOPE->nf.mVector.push_back(newM);
            }
        }
    }


}

static void lSpheres(string sphereList)
{
    /*This function loads a list of spheres given in the string sphereList
        The format is:
            x y z a m
        where
            x, y, z = sphere position (required)
            a = sphere radius (required)
            m = material ID (optional) */

    std::stringstream ss(sphereList);

    while(!ss.eof())
    {
        //create a new sphere
        sphere newS;

        //get the sphere data
        ss>>newS.p[0];
        ss>>newS.p[1];
        ss>>newS.p[2];
        ss>>newS.a;

        if(ss.peek() != '\n')
            ss>>newS.iMaterial;

        //add the new sphere to the sphere vector
        SCOPE->nf.sVector.push_back(newS);

        //ignore the rest of the line
        ss.ignore(1000, '\n');

        //check out the next element (this should set the EOF error flag)
        ss.peek();
    }
}

void lSpheres(rts::arglist args)
{
    //if a sphere is specified at the command line
    if(args("spheres"))
    {
        rts::argument sphere_arg = args["spheres"];
		int nArgs = sphere_arg.nargs();

        //compute the number of spheres specified
        unsigned int nS;
        if(nArgs <= 5)
            nS = 1;
        else
        {
            //if the number of parameters is divisible by 4, compute the number of spheres
            if(nArgs % 5 == 0)
                nS = nArgs / 5;
            else
            {
                cout<<"BIMSIM Error: Invalid number of sphere parameters."<<endl;
                exit(1);
            }
        }

        stringstream ss;

        //for each sphere
        for(unsigned int s=0; s<nS; s++)
        {
            //compute the number of sphere parameters
            unsigned int nP;
            if(nS == 1) nP = nArgs;
            else nP = 5;

            //store each parameter as a string
            for(unsigned int i=0; i<nP; i++)
            {
                ss<<sphere_arg.as_float(s*5 + i)<<" ";
            }
            ss<<endl;
        }



        //convert the string to a sphere list
        lSpheres(ss.str());
    }

    //if a files are specified
    if(args("sphere-file"))
    {
		rts::argument sfiles = args["sphere-file"];
		int nFiles = sfiles.nargs();

        //load each file
        for(unsigned int iS=0; iS<nFiles; iS++)
        {
            //load the file into a string
            std::ifstream ifs(sfiles.as_text(iS).c_str());

            if(!ifs)
            {
                cout<<"Error loading sphere file."<<endl;
                exit(1);
            }

            std::string instr((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

            //load the list of spheres from a string
            lSpheres(instr);
        }
    }

    //make sure the appropriate materials are loaded
    unsigned int nS = SCOPE->nf.sVector.size();

    //for each sphere
    for(unsigned int s = 0; s<nS; s++)
    {
      //  std::cout<<"\t\t SCOPE->nf.mVector.size(): "<<SCOPE->nf.mVector.size()<<std::endl;
        //make sure the corresponding material exists
        if(SCOPE->nf.sVector[s].iMaterial + 1 > SCOPE->nf.mVector.size())
        {
            //otherwise output an error
            cout<<"BIMSIM Error - A material is not loaded for sphere "<<s+1<<"."<<endl;
            cout<<"Material requested: "<<SCOPE->nf.sVector[s].iMaterial + 1<<endl;
            cout<<"Number of materials: "<<SCOPE->nf.mVector.size()<<endl;
            exit(1);
        }
    }
}

static void lOptics(rts::arglist &args)
{
    if(args("objective"))
    {
        
        if(args["objective"].nargs() == 1)
		{
			SCOPE->objective[0] = 0.0;
            SCOPE->objective[1] = args["objective"].as_float();
		}
        else
        {
            SCOPE->objective[0] = args["objective"].as_float(0);
            SCOPE->objective[1] = args["objective"].as_float(1);
        }
    }
}

static void lImagePlane(rts::arglist args)
{
	bsPoint pMin(DEFAULT_PLANE_MIN_X, DEFAULT_PLANE_MIN_Y, DEFAULT_PLANE_MIN_Z);
	bsPoint pMax(DEFAULT_PLANE_MAX_X, DEFAULT_PLANE_MAX_Y, DEFAULT_PLANE_MAX_Z);
	bsVector normal(DEFAULT_PLANE_NORM_X, DEFAULT_PLANE_NORM_Y, DEFAULT_PLANE_NORM_Z);

    //default plane size in microns
    ptype s = DEFAULT_PLANE_SIZE;
    ptype pos = DEFAULT_PLANE_POSITION;

	if(args("xy"))
	{


		
		if(args["xy"].nargs() >= 1)
			s = args["xy"].as_float(0);
		if(args["xy"].nargs() == 2)
			pos = args["xy"].as_float(1);

        cout<<"\t\t pos: "<<pos<<endl;

		//calculate the plane corners and normal based on the size and position
		pMin = bsPoint(-s/2, -s/2, pos);
		pMax = bsPoint(s/2, s/2, pos);
		normal = bsVector(0, 0, 1);
	}
	else if(args("xz"))
	{
		//default plane size in microns
		ptype size = DEFAULT_PLANE_SIZE;

		if(args["xz"].nargs() >= 1)
			size = args["xz"].as_float(0);
		if(args["xz"].nargs() >= 2)
			pos = args["xz"].as_float(1);

		//calculate the plane corners and normal based on the size and position
		pMin = bsPoint(-size/2, pos, -size/2);
		pMax = bsPoint(size/2, pos, size/2);
		normal = bsVector(0, -1, 0);
	}
	else if(args("yz"))
	{
		//default plane size in microns
		ptype size = DEFAULT_PLANE_SIZE;

		if(args["yz"].nargs() >= 1)
			size = args["yz"].as_float(0);
		if(args["yz"].nargs() >= 2)
			pos = args["yz"].as_float(1);

		//calculate the plane corners and normal based on the size and position
		pMin = bsPoint(pos, -size/2, -size/2);
		pMax = bsPoint(pos, size/2, size/2);
		normal = bsVector(1, 0, 0);
	}
	//set the default values for the slice position and orientation
	else{


		pMin = bsPoint(args["plane-lower-left"].as_float(0), args["plane-lower-left"].as_float(1), args["plane-lower-left"].as_float(2));
		pMax = bsPoint(args["plane-upper-right"].as_float(0), args["plane-upper-right"].as_float(1), args["plane-upper-right"].as_float(2));
		normal = bsVector(args["plane-normal"].as_float(0), args["plane-normal"].as_float(1), args["plane-normal"].as_float(2));

	}
	
	SCOPE->setPos(pMin, pMax, normal);

    SCOPE->setZdPos(pos);
    //cout<<"\t\t lambda: "<<SCOPE->nf.lambda<<endl;
    float dx = s/args["resolution"].as_float();
    float lambda = SCOPE->nf.lambda;

    int ss = 1;
    if (dx > lambda/2.0)
        ss = (int) ceil(2*s/(args["resolution"].as_float()*lambda));

    cout<<"\t\t ss: "<<ss<<endl;

	//resolution
//	SCOPE->setRes(args["resolution"].as_float(),
//				  args["resolution"].as_float(),
//				  args["padding"].as_float(),
//				  args["supersample"].as_float());

    //resolution
    SCOPE->setRes(args["resolution"].as_float(),
                  args["resolution"].as_float(),
                  args["padding"].as_float(),
                  ss);






	SCOPE->setNearfield();
}

static void lNearfield(rts::arglist args)
{
    //test to see if we are running a vector field simulation
    bool vectorField = false;
    if(args("vector"))
        vectorField = true;
    SCOPE->scalarSim = !vectorField;

	//test to see if we are simulating a plane wave
	bool planeWave = DEFAULT_PLANEWAVE;
	if(args("plane-wave"))
		planeWave = !planeWave;
	SCOPE->nf.planeWave = planeWave;

	//get the incident field amplitude
	SCOPE->nf.A = args["amplitude"].as_float();

	//get the condenser parameters

    if(args["condenser"].nargs() == 1)
	{
        SCOPE->nf.condenser[0] = 0;
		SCOPE->nf.condenser[1] = args["condenser"].as_float(0);
	}
    else
    {
        SCOPE->nf.condenser[0] = args["condenser"].as_float(0);
        SCOPE->nf.condenser[1] = args["condenser"].as_float(1);
    }


    //get the extended source parameters

    SCOPE->nf.extSource[0] = args["points-rings"].as_int(0);
    SCOPE->nf.extSource[1] = args["points-rings"].as_int(1);

   // cout<<"\t\t\t nf.extSource = "<<SCOPE->nf.extSource[0]<<" , "<<SCOPE->nf.extSource[1]<<endl;

	//get the focal rtsPoint position
    SCOPE->nf.focus[0] = args["focus"].as_float(0);
    SCOPE->nf.focus[1] = args["focus"].as_float(1);
    SCOPE->nf.focus[2] = args["focus"].as_float(2);

	//get the incident light direction (k-vector)
	bsVector spherical(1, 0, 0);

    //if a k-vector is specified
    if(args("k"))
    {
        
        spherical[1] = args["k"].as_float(0);
        spherical[2] = args["k"].as_float(1);
    }
	SCOPE->nf.k = spherical.sph2cart();


    //incident field order
    SCOPE->nf.m = args["field-order"].as_int();

    //number of Monte-Carlo samples
    SCOPE->nf.nWaves = args["samples"].as_int();

	//random number seed for Monte-Carlo samples
	if(args("seed"))
		srand(args["seed"].as_int());
}

static void lOutputParams(rts::arglist args)
{
    //append simulation results to previous binary files
    gFileOut.append = DEFAULT_APPEND;
    if(args("append")){

        gFileOut.append = true;
    //    std::cout<<"\t\t append is on"<<std::endl;
    }



	//image parameters
	//component of the field to be saved
	std::string fieldStr;
    fieldStr = args["output-type"].as_text();

    if(fieldStr == "magnitude")
        gFileOut.field = fileoutStruct::fieldMag;
    else if(fieldStr == "intensity")
        gFileOut.field = fileoutStruct::fieldIntensity;
    else if(fieldStr == "incident")
        gFileOut.field = fileoutStruct::fieldIncident;
    else if(fieldStr == "polarization")
        gFileOut.field = fileoutStruct::fieldPolar;
    else if(fieldStr == "imaginary")
        gFileOut.field = fileoutStruct::fieldImag;
    else if(fieldStr == "real")
        gFileOut.field = fileoutStruct::fieldReal;
    else if(fieldStr == "angular-spectrum")
        gFileOut.field = fileoutStruct::fieldAngularSpectrum;


	//image file names
	gFileOut.intFile = args["intensity"].as_text();
    gFileOut.incFile = args["incident"].as_text();

	gFileOut.absFile = args["absorbance"].as_text();

    gFileOut.absSpecFile = args["absorbance-spec"].as_text();

	if(args("transmittance"))
		gFileOut.transFile = args["transmittance"].as_text();
	gFileOut.nearFile = args["near-field"].as_text();
	gFileOut.farFile = args["far-field"].as_text();

	//colormap
	std::string cmapStr;
    cmapStr = args["colormap"].as_text();
    if(cmapStr == "brewer")
        gFileOut.colormap = stim::cmBrewer;
    else if(cmapStr == "gray")
        gFileOut.colormap = stim::cmGrayscale;
    else
        cout<<"color-map value not recognized (using default): "<<cmapStr<<endl;
}

void LoadParameters(rts::arglist &args)
{
    lFlags(args);
    lWavelength(args);
	lMaterials(args);
	lSpheres(args);
	lOptics(args);
	lImagePlane(args);
	lNearfield(args);
	lOutputParams(args);


	//if an extended source will be used
    if(args("extended-source"))
    {
        //load the point sources
        std::string filename = args["extended-source"].as_text();
        SCOPE->LoadExtendedSource(filename);

    }else if(args("interpolate")) //if an extended source will be used with interpolation
    {

        SCOPE->interpolate = args["interpolate"].as_int();
//        if(SCOPE->interpolate==1)
//            //generate one point source at [0 0 a]
//            SCOPE->generatePointSources(1,1,(scope->nf.sVector[1]).a);
    }


}

static void OutputOptions()
{
	cout<<SCOPE->toStr();

    //cout<<"# of source points: "<<SCOPE->focalPoints.size()<<endl;

}

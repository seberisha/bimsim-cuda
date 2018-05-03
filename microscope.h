#ifndef DETECTORSTRUCT_H
#define DETECTORSTRUCT_H

#include "fieldslice.h"
#include "nearfield.h"

#include <iostream>
using namespace std;

struct sourcePoint
{
    bsPoint f;
    ptype A;
};

struct sourcePointES
{
    bsPoint f;
    ptype A;
    ptype rotationAngle;
};

struct microscopeStruct
{
    //int interpolate=0;
    int interpolate;
    int padding;

    int ss;

    //float scale = 1.0;
    float scale;

    //field slice at the detector
    fieldslice Ufd;
    fieldslice Ud;
    bool scalarSim;

    //detector image
    scalarslice* D;		//sample

    scalarslice* D_ps1;		//detector image of the total field for one point source
    scalarslice* Di_ps1;    //detector image of the focues field for one point source

    ptype sumD; //sum of elements of D
    scalarslice* Di;		//incident field image
    ptype sumDi;    //sum of elements of Di

    ptype objective[2];

    //image position and orientation in world space
    rts::quad<ptype, 3> pos;
    float z_d;

    vector<sourcePoint> focalPoints;
    vector<sourcePointES> focalPointsES;

    microscopeStruct();

    nearfieldStruct nf;


    void setZdPos(float d);

    //output results for each point source
    void saveDetectorForEachPointSource();



    //save the total radius of the final detector image to
    //texture
    void SaveLineToTexture(int radiusLength);

    void InterpLineAtDetector();


    void SimulateExtendedSourceInterpLineAtDetector(int *pointsRings);


    //interpolate
    void InterpLineForAllPointSources(int numSamples, int numRings);

    //simulate Gaussian extended source without interpolation
    void SimulateExtendedGausssianSourceNoInterp(int *pointsRings);


    void SaveDetectorToTexture();

    void init();
    void destroy();

    //functions for dealing with extended sources
    void SimulateExtendedSource();

    void SimulateExtendedSourceWithInterpolation(float sphereRadius, int *pointsRings);

    //simulate extended source by interpolating fields at the detector
    void SimulateExtendedSourceInterpAtDetector(int *pointsRings);
    //interpolate the total and focused fields at the detector
    void InterpAtDetector();

    void saveDetectorForOnePs();

    void LoadExtendedSource(std::string filename);

    void generatePointSources(int numPoints, float radius, float amplitude);
    void generatePointSourcesNoInterp(float  spacing, int radius, float amplitude);

    void SimulateScattering();

    void SimulateInterpolation();

    void SimulateImaging();
    void Simulate();
    void SimulateAndSaveTexture();
    void SimulateWithInterpolation();

    void saveToTexture();

    void setNearfield();
    void setRes(int x_res, int y_res, int pad, int supersampling);
    void setPos(bsPoint pMin, bsPoint pMax, bsVector normal);
	
    void applyBackpropagation();
    void applyBandpass();
    void getFarField();
    void integrateDetector();
    void clearDetector();

    //compute detector measurements
    scalarslice getAbsorbance();

    ptype getAbsorbanceSpectrum();

    scalarslice getTransmittance();
    scalarslice getIntensity();
    scalarslice getIncidentFieldImage();
    string toStr();



};

#endif

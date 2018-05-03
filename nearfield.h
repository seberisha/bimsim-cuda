#ifndef NEARFIELD_H
#define NEARFIELD_H

//#include "defaults.h"
#include "fieldslice.h"
#include "montecarlo.h"
#include "rts/optics/material.h"
#include "sphere.h"
#include <vector>

#define EPSILON_FLOAT   0.000001f

//This structure stores values relevant to creating the near field
struct nearfieldStruct
{
	//incident wavelength
	ptype lambda;

	//condenser numerical aperature (internal and external)
	ptype condenser[2];
    int extSource[2];

	//amplitude of the incident field
	ptype A;

	//incident field polarization;
	bsVector E;

	//position of the focus in 3D space
	bsVector k;		//cartesian coordinates, normalized
	bsPoint focus;
    ptype rotationAngle;

	//slice position and orientation in world space
	rts::quad<ptype, 3> pos;

	//slices for the focused field
	fieldslice Uf;
	ptype d_min, d_max;

	//	and total field: Uf + sum(Us)
	fieldslice U;

	//incident field order
	int m;

	//flag for a vector simulation
	bool scalarSim;

	//flag for a plane wave
	bool planeWave;

	//flag for using a LUT
	bool lut_uf;
	bool lut_us;

	//timings
	float t_Uf;
	float t_Us;



	//---------Scatterers------------

	//angular resolution of scatterers
	//		number of divisions in [0 pi]
	unsigned int angRes;

	//list of materials
	std::vector< rts::material<ptype> > mVector;

	//list of scatterers
	std::vector<sphere> sVector;

	nearfieldStruct();

    void saveFieldsToTexture();
    void Interpolate();

	void init();
	void destroy();

	void Simulate();

	void calcSpheres();

	//plane waves for Monte-Carlo sampling
	std::vector<bsVector> inWaves;
	int nWaves;
	void calcWaves();

	std::string toStr();

	void setRes(int x_res, int y_res);

	void setPos(bsPoint pMin, bsPoint pMax, bsVector normal);

	//this function re-computes the focused field
	void calcUf();

	//scalar functions
	void scalarUf();
	void scalarUfLut();

	//vector functions
	void vectorUf();
	void vectorUfLut();

	void calcBesselLut(ptype* j, ptype d_min, ptype d_max, int dR);

	//compute the field scattered by all of the materials
	void calcUs();
	void scalarUs();
	void scalarUpLut();


	//add the incident field to the sum of scattered fields
	void sumUf();

};



#endif

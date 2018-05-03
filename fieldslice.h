#ifndef FIELDSLICE_H
#define FIELDSLICE_H

#include "dataTypes.h"

#include <string>
#include <sstream>
#include "rts/math/quad.h"

#include "scalarslice.h"


struct fieldslice
{
	//images representing the slice for each world-space polarization direction
	bsComplex* x_hat;
	bsComplex* y_hat;
	bsComplex* z_hat;



	//slice resolution
	unsigned int R[2];

	//if the slice is a scalar
	bool scalarField;



	fieldslice(unsigned int x_res, unsigned int y_res);

	fieldslice();

	~fieldslice();

	//copy constructor
	fieldslice(const fieldslice& rhs);

	//void setPos(bsPoint pMin, bsPoint pMax, bsVector N);

	scalarslice Mag();
	scalarslice Real();
	scalarslice Imag();
	scalarslice Intensity();
    void IntegrateAndResample(scalarslice* detector, unsigned int supersample, float scale);

	void ScaleField(ptype v);

	//convert the field slice to the angular spectrum
	void toAngularSpectrum();
	void fromAngularSpectrum();

	//crop a region from the field
	fieldslice crop(int u, int v, int su, int sv);
	fieldslice& operator=(const fieldslice& rhs);

	void init_gpu();
	void kill_gpu();
	void clear_gpu();

	std::string toStr();

    //save the detector images for one point source so we can use them to interpolate for the rest
    //of the point sources
    void ResampleAndSave(scalarslice* detector, unsigned int supersample);


};

#endif

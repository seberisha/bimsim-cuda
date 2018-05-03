#include "montecarlo.h"
#include "rts/math/quaternion.h"
#include "rts/math/matrix.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

void mcSampleNA(bsVector* samples, int N, bsVector k, ptype NAin, ptype NAout)
{
    /*Create Monte-Carlo samples of a cassegrain objective by performing uniform sampling
        of a sphere and projecting these samples onto an inscribed sphere.

        samples = rtsPointer to sample vectors specified as normalized cartesian coordinates
        N       = number of samples
		kSph	= incident light direction in spherical coordinates
        NAin    = internal obscuration NA
        NAout   = outer cassegrain NA
    */

	//get the axis of rotation for transforming (0, 0, 1) to k
	//k = -k;
	ptype cos_angle = k.dot(bsVector(0, 0, 1));

//    if(verbose)
//    {
//        cout<<"monteCarlo====K Vector:"<<"["<<k[0]<<" "<<k[1]<<" "<<k[2]<<"]"<<endl;

//    }

	rts::matrix<ptype, 3> rotation;
	if(cos_angle != 1.0)
	{
		bsVector axis = bsVector(0, 0, 1).cross(k).norm();

		ptype angle = acos(cos_angle);
		rts::quaternion<ptype> quat;
		quat.CreateRotation(angle, axis);
		rotation = quat.toMatrix3();
	}

    //find the phi values associated with the cassegrain ring
    ptype inPhi = asin(NAin);
    ptype outPhi = asin(NAout);

    //calculate the z-values associated with these angles
    ptype inZ = cos(inPhi);
    ptype outZ = cos(outPhi);

    ptype rangeZ = inZ - outZ;

    //draw a distribution of random phi, z values
    ptype z, phi, theta;
    for(int i=0; i<N; i++)
    {
        z = ((ptype)rand() / (ptype)RAND_MAX) * rangeZ + outZ;
        theta = ((ptype)rand() / (ptype)RAND_MAX) * 2 * (ptype)PI;

        //calculate theta
        phi = acos(z);

        //compute and store cartesian coordinates
		bsVector spherical(1, theta, phi);
		bsVector cart = spherical.sph2cart();
        samples[i] = rotation * cart;
    }




}

#ifndef RTS_SCALAR_SLICE
#define RTS_SCALAR_SLICE

#include "dataTypes.h"
#include "stim/visualization/colormap.h"

struct scalarslice
{
	//gpu pointer to the scalar slice
	ptype* S;

    //gpu pointer to the radius of the scalar slice
    ptype *line;

	//resolution of the slice
	int R[2];

    ptype rotationAngle;

    scalarslice();
	scalarslice(int x, int y);
    scalarslice(int x, int y, int psSamples);
    scalarslice(ptype s);
	~scalarslice();
	void clear();

    void clearLine(int radiusLength);

	void toImage(std::string filename, ptype vmin, ptype vmax, stim::colormapType cmap = stim::cmBrewer);
	void toImage(std::string filename, bool positive = true, stim::colormapType cmap = stim::cmBrewer);
	void toEnvi(std::string filename, ptype wavelength = 0, bool append = false);

	//assignment operator
	scalarslice & operator= (const scalarslice & rhs);

    void SaveDetectorToTexture(scalarslice *D_ps1, scalarslice*Di_ps1);

};



#endif

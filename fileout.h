#ifndef FILE_OUTPUT_H
#define FILE_OUTPUT_H

#include <string>
#include "dataTypes.h"

#include "stim/visualization/colormap.h"
#include "fieldslice.h"
#include "nearfield.h"
#include "microscope.h"
#include "rts/cuda/error.h"

struct fileoutStruct{

	//output file names
	std::string nearFile;		//near field filename
	std::string farFile;		//far field filename
	std::string intFile;		//detector intensity filename
    	std::string incFile;		//detector incident filename

	std::string absFile;        //detector absorbance filename

    	std::string absSpecFile;        //detector absorbance spectrum filename

	std::string transFile;      //detector transmission filename

	bool append;                //append simulation results to existing binary files

	//output type
    	enum field_type {fieldMag, fieldIntensity, fieldIncident, fieldAbsorbance, fieldAbsorbanceSpectrum,fieldPolar, fieldImag, fieldReal, fieldAngularSpectrum};
	enum image_source {imageNearField, imageObjective, imageDetector, imageExtendedSource};

	field_type field;

	//flag for output in wavenumber
	bool wavenumber;

	//image_source source;

	//color map info
	stim::colormapType colormap;
	ptype colorMax;

	void Save(microscopeStruct* scope);
	void Simulate(microscopeStruct* scope);

	private:
	bool is_binary(std::string filename);
	void saveNearField(nearfieldStruct* nf);
	void saveFarField(microscopeStruct* scope);
	void saveDetector(microscopeStruct* scope);

};



#endif

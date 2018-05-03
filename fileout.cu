#include "fileout.h"

void fileoutStruct::saveNearField(nearfieldStruct* nf)
{
    if(nearFile == "") return;

    if(field == fieldReal)
    {
        scalarslice S = nf->U.Real();
        S.toImage(nearFile, false, colormap);
    }
    if(field == fieldImag)
    {
        scalarslice S = nf->U.Imag();
        S.toImage(nearFile, false, colormap);
    }
    if(field == fieldMag)
    {
        scalarslice S = nf->U.Mag();
        S.toImage(nearFile, true, colormap);
    }
}


void fileoutStruct::saveFarField(microscopeStruct* scope)
{
    if(farFile == "") return;

    if(field == fieldReal)
    {
        scalarslice S = scope->Ud.Real();
        S.toImage(farFile, false, colormap);
    }
    if(field == fieldImag)
    {
        scalarslice S = scope->Ud.Imag();
        S.toImage(farFile, false, colormap);
    }
    if(field == fieldMag)
    {
        scalarslice S = scope->Ud.Mag();
        S.toImage(farFile, true, colormap);
    }

}

void fileoutStruct::saveDetector(microscopeStruct* scope)
{

    //intensity
    if(intFile != "")
    {
        cout<<"\t\t saving: "<<intFile<<endl;
        scalarslice I = scope->getIntensity();

        if(is_binary(intFile))
        {
            if(wavenumber)
                I.toEnvi(intFile, 10000.0f/scope->nf.lambda, append);
            else
                I.toEnvi(intFile, scope->nf.lambda, append);
        }
        else
            //intFile = intFile + ".bmp";
            I.toImage(intFile);
    }

    //incident field image

    //intensity
    if(incFile != "")
    {
        cout<<"\t\t saving: "<<incFile<<endl;
        scalarslice I = scope->getIncidentFieldImage();
        //  std::cout<<"\t\t=====> incident file: "<<incFile<< std::endl;

        if(is_binary(incFile))
        {
            if(wavenumber)
                I.toEnvi(incFile, 10000.0f/scope->nf.lambda, append);
            else
                I.toEnvi(incFile, scope->nf.lambda, append);
        }
        else
            //incFile = incFile + ".bmp";
            I.toImage(incFile);
    }


    //std::cout<<"\t\t******=====> writing absorbance file: "<<std::endl;
    //absorbance
    if(absFile != "")
    {
        cout<<"\t\t saving: "<<absFile<<endl;
        scalarslice I = scope->getAbsorbance();

        if(is_binary(absFile))
        {
            if(wavenumber)
                I.toEnvi(absFile, 10000.0f/scope->nf.lambda, append);
            else
                I.toEnvi(absFile, scope->nf.lambda, append);
        }
        else
        {
            //cout<<"toImage(absFile)"<<endl;
            //absFile = absFile + ".bmp";
            I.toImage(absFile);
        }
    }

    //std::cout<<"\t\t******=====> writing absorbance spectrum file: "<<std::endl;

    //absorbance
    if(absSpecFile != "")
    {
        cout<<"\t\t saving: "<<absSpecFile<<endl;
        ptype a = scope->getAbsorbanceSpectrum();

        // std::cout<<"\t\t=====> a is: "<<a<<std::endl;



        if(is_binary(absSpecFile))
        {
            ofstream myfile;
            myfile.open (absSpecFile.c_str(), ios::app);

            if(wavenumber)
                myfile << 10000.0f/scope->nf.lambda <<","<<a<<endl;
            else
                myfile << scope->nf.lambda <<","<<a<<endl;
            myfile.close();


        }
        //  else
        //  I.toImage(absSpecFile);
    }


    //transmittance
    if(transFile != "")
    {
        cout<<"\t\t saving: "<<transFile<<endl;
        scalarslice I = scope->getTransmittance();

        if(is_binary(transFile))
        {
            if(wavenumber)
                I.toEnvi(transFile, 10000.0f/scope->nf.lambda, append);
            else
                I.toEnvi(transFile, scope->nf.lambda, append);
        }
        else
            I.toImage(transFile);
    }

}

bool fileoutStruct::is_binary(std::string filename)
{
    //this function guesses if a file name is binary or a standard image
    //  do this by just testing extensions

    //get the extension
    size_t i = filename.find_last_of('.');

    //if there is no extension, return true
    if( i == std::string::npos )
        return true;

    //otherwise grab the extension
    std::string ext = filename.substr(i+1);
    if(ext == "bmp" ||
            ext == "jpg" ||
            ext == "tif" ||
            ext == "gif" ||
            ext == "png")
        return false;
    else
        return true;
}



void fileoutStruct::Save(microscopeStruct* scope)
{


    //save images of the fields in the microscope

    //printf("\n\t\t\t Save\n");


    //if the user specifies an extended source


    //if the user specifies an extended source
    if(scope->focalPoints.size() > 0)
    {
        std::cout<<"\t\t ====>extended source simulation with image"<<std::endl;
        //simulate the extended source and output the detector image
        scope->SimulateExtendedSource();

        // saveNearField(&scope->nf);
        //saveFarField(scope);

        //save the detector images
        saveDetector(scope);

        //simulate scattering for the last point (so that you have a near field image)
        //scope->SimulateScattering();
        //saveNearField(&scope->nf);

    }else
    {
        //interpolate the nearfield
        if (scope->interpolate == 3)
        {
            scope->SimulateExtendedSourceWithInterpolation(0, scope->nf.extSource);

            saveDetector(scope);
        }
        else if (scope->interpolate == 1)
        {
            //interpolate intensity at the detector
            scope->SimulateExtendedSourceInterpAtDetector(scope->nf.extSource);
            saveDetector(scope);
        }else if (scope->interpolate == 0)
        {
            std::cout<<"\t\t ====>scope->interpolate == 0 "<<std::endl;
            scope->SimulateExtendedGausssianSourceNoInterp(scope->nf.extSource);

            //save the detector images
            saveDetector(scope);
        }else if (scope->interpolate == 2)
        {
            // interpolate one line at the detector
            std::cout<<"\t\t ====>scope->interpolate == 2 "<<std::endl;
            scope->SimulateExtendedSourceInterpLineAtDetector(scope->nf.extSource);
            saveDetector(scope);
        }else
        {
            //std::cout<<"\t\t ====>single source simulation -- SimulateScattering"<<std::endl;
            //run the near-field simulation
            //scope->SimulateScattering();

            //std::cout<<"\t\t========================= SimulateImaging "<<std::endl;
            //output the near field image
            //saveNearField(&scope->nf);

            //run the far-field simulation
            //scope->SimulateImaging();
	        std::cout<<"\t\t single source simulation -- original version"<<endl;
	        scope->Simulate();
            saveNearField(&scope->nf);
            saveFarField(scope);
            std::cout<<"\t\t saveDetector"<<std::endl;
            //save the detector images
            saveDetector(scope);

        }

    }


}




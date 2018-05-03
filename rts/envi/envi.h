#ifndef RTS_ENVI_H
#define RTS_ENVI_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "rts/envi/envi_header.h"

/*	This is a class for reading and writing ENVI binary files. This class will be updated as needed.

What this class CAN currently do:
	*)  Write band images to a BSQ file.
	*)  Append band images to a BSQ file.
	*)  Interpolate and query bands in a float32 BSQ file
*/
class EnviFile
{
	EnviHeader header;
	FILE* data;
	std::string mode;

	//a locked file has alread had data written..certain things can't be changed
	//		accessor functions deal with these issues
	bool locked;

	void init()
	{
		locked = false;
		data = NULL;
	}

	void readonly()
	{
		if(mode == "r")
		{
			std::cout<<"ENVI Error: changes cannot be made to a read-only file."<<std::endl;
			exit(1);
		}
	}

	//this function determines if the current software is capable of the requested change
	bool caps(EnviHeader new_header)
	{
		readonly();		//if the file is read-only, throw an exception

		//we currently only support BSQ
		if(new_header.interleave == Envi::BIP || new_header.interleave == Envi::BIL)
		{
			std::cout<<"This software only supports BSQ.";
			return false;
		}

		//if the number of bands has changed
		if(header.bands != new_header.bands)
		{
			//the new header must be a BSQ file
			if(new_header.interleave != Envi::BSQ)
			{
				std::cout<<"ENVI Error: can only add bands to a BSQ file."<<std::endl;
				return false;
			}
		}

		//disallow anything that can't be done when the file is locked
		if(locked)
		{
			if(header.lines != new_header.lines)
			{
				std::cout<<"ENVI Error: cannot change the number of lines from "<<header.lines<<" to "<<new_header.lines<<" in a locked BSQ file."<<std::endl;
				return false;
			}
			if(header.samples != new_header.samples)
			{
				std::cout<<"ENVI Error: cannot change the number of samples in a locked BSQ file."<<std::endl;
				return false;
			}
			if(header.data_type != new_header.data_type)
			{
				std::cout<<"ENVI Error: cannot change the data type of a locked file."<<std::endl;
				return false;
			}
			if(header.byte_order != new_header.byte_order)
			{
				std::cout<<"ENVI Error: cannot change the byte order of a locked file."<<std::endl;
				return false;
			}
		}

		return true;

	}

	EnviHeader load_header(std::string header_file)
	{
		EnviHeader new_header;
		new_header.load(header_file);
		return new_header;
	}

	void set_header(EnviHeader new_header)
	{
		if(caps(new_header))
			header = new_header;
		else
			exit(1);
	}

	bool test_exist(std::string filename)
	{
        //tests for the existance of the specified binary file
        FILE* f;
        long sz = 0;
        f = fopen(filename.c_str(), "rb");
        if(f == NULL)
        {
            //std::cout<<"ENVI Error: error testing file existance."<<std::endl;
            //exit(1);
			return false;
        }

        fseek(f, 0, SEEK_END);
        sz = ftell(f);
        fclose(f);

        if(sz>0)
            return true;
        else
            return false;

    }

    void get_surrounding_bands(double wavelength, unsigned int &low, unsigned int &high)
    {
        int i;
        int nBands = header.wavelength.size();
        low = high = 0;
        //std::cout<<"Wavelength to search: "<<wavelength<<"     Number of Bands: "<<nBands<<std::endl;
        for(i=0; i<nBands; i++)
        {
            if(header.wavelength[i] < wavelength)
            {
                if(header.wavelength[low] > wavelength)
                    low = i;
                if(header.wavelength[i] > header.wavelength[low])
                    low = i;
            }

            if(header.wavelength[i] > wavelength)
            {
                if(header.wavelength[high] < wavelength)
                    high = i;
                if(header.wavelength[i] < header.wavelength[high])
                    high = i;
            }
            //std::cout<<"Low: "<<low<<"     High: "<<high<<std::endl;
        }
        //exit(1);
    }

    void interpolate_band(void* result, void* A, void* B, double a)
    {
        //interpolate between two bands independent of data type
        //  CURRENTLY ONLY HANDLES FLOAT
        int X = header.samples;
        int Y = header.lines;
        int N = X * Y;

        if(header.data_type != Envi::float32)
        {
            std::cout<<"ENVI Error: this software only handles interpolation of 32-bit floats."<<std::endl;
            exit(1);
        }

		float v0, v1;
        for(int n=0; n<N; n++)
        {
			v0 = ((float*)A)[n];
			v1 = ((float*)B)[n];
            float r = v0 * a + v1 * (1.0 - a);
            //((float*)result)[n] = ((float*)A)[n] * a + ((float*)B)[n] * (1.0 - a);
            ((float*)result)[n] = r;
        }

    }

public:

	EnviFile()
	{
        mode = "w";
		init();
	}

	EnviFile(std::string filename, std::string file_mode = "r")
	{
		init();
		open(filename, filename + ".hdr", file_mode);
	}

	void open(std::string file_name, std::string header_name, std::string file_mode = "r")
	{

        //load the header file
        //  if the file doesn't exist, that is okay unless it is read-only

        //if the file is being written, create a new header
        if(file_mode == "w")
            header.name = header_name;
        else if(!header.load(header_name))
        {
            if(file_mode == "r")
            {
                std::cout<<"ENVI Error: header file not found: "<<header_name<<std::endl;
                exit(1);
            }
            //otherwise use the default header and save the name
            header.name = header_name;
        }

		mode = file_mode;

		//lock the file if it is read-only
		if(mode == "r")
            locked = true;
		else if(mode == "a" || mode == "r+" || mode == "a+")
		{
            //if the file can be edited, lock it if data has already been written
            if(test_exist(file_name))
                locked = true;
		}


		//open the data file
		std::string tmpmode = mode+"b";
		//std::cout<<"Mode: "<<tmpmode<<std::endl;
		data = fopen(file_name.c_str(), tmpmode.c_str());
		if(data == NULL)
		{
			std::cout<<"ENVI Error: unable to open binary file: "<<file_name<<std::endl;
			exit(1);
		}


	}

	void setDescription(std::string desc)
	{
		readonly();		//if the file is read-only, throw an exception
		header.description = desc;
	}

	void setZPlotTitles(std::string spectrum, std::string value)
	{
		readonly();		//if the file is read-only, throw an exception
		header.z_plot_titles[0] = spectrum;
		header.z_plot_titles[1] = value;
	}

	void setPixelSize(double sx, double sy)
	{
		readonly();		//if the file is read-only, throw an exception
		header.pixel_size[0] = sx;
		header.pixel_size[1] = sy;
	}

	void setWavelengthUnits(std::string units)
	{
		readonly();		//if the file is read-only, throw an exception
		header.wavelength_units = units;
	}

	void setCoordinates(double x, double y)
	{
		readonly();		//if the file is read-only, throw an exception
		header.x_start = x;
		header.y_start = y;
	}

	//FUNCTIONS THAT MAY BE DISALLOWED IN A LOCKED FILE
	void setSamples(unsigned int samples)
	{
		EnviHeader newHeader = header;
		newHeader.samples = samples;
		set_header(newHeader);
	}

	void setLines(unsigned int lines)
	{
		EnviHeader newHeader = header;
		newHeader.lines = lines;
		set_header(newHeader);
	}

	unsigned int lines() {return header.lines;}
	unsigned int samples() {return header.samples;}
	unsigned int bands() {return header.bands;}
	double wavelength(unsigned int b) {return header.wavelength[b];}

	void setBands(unsigned int bands)
	{
		EnviHeader newHeader = header;
		newHeader.bands = bands;
		set_header(newHeader);
	}

	unsigned int getElementSize()
	{
        //return the size of each element in bytes
        return header.valsize();
    }
    unsigned int getBandSize()
    {
        //returns the size of a band in bytes
        return header.valsize() * header.lines * header.samples;
    }

    Envi::dataType getDataType()
    {
        return header.data_type;
    }

	//DATA MANIPULATION
	void addBand(void* band, unsigned int samples, unsigned int lines, double wavelength, std::string band_name ="")
	{
		//add a band to the file

		EnviHeader newHeader = header;
		newHeader.bands++;
		newHeader.samples = samples;
		newHeader.lines = lines;
		newHeader.wavelength.push_back(wavelength);
		if(band_name != "")
			newHeader.band_names.push_back(band_name);
		set_header(newHeader);

		//copy the data to the file
		fwrite(band, header.valsize(), samples * lines, data);

		//since data was written, the file must be locked
		locked = true;
	}

	//DATA RETRIEVAL
	void getBand(void* ptr, unsigned int n)
	{
        //reads the n'th band and writes it to the given pointer
        //  memory has been assumed to be allocated

        if(n > header.wavelength.size())
        {
            std::cout<<"ENVI Error: Invalid band number."<<std::endl;
            return;
        }


        if(header.interleave == Envi::BSQ)
        {

            //compute the position of the desired band (relative to the beginning of the file)
            int bandpos = header.header_offset + getBandSize() * n;

            //seek to that position
            int seek_result = fseek(data, bandpos, SEEK_SET);
            if(seek_result)
            {
                std::cout<<"ENVI Error: Error seeking through data file."<<std::endl;
                return;
            }

            //read the band data from the file
            size_t n_read = fread(ptr, getElementSize(), header.samples * header.lines, data);
            if(n_read != header.samples * header.lines)
            {
                std::cout<<"ENVI Error -- Error reading band data: "<<n_read<<"/"<<samples() * lines()<<" pixels read"<<std::endl;
                exit(1);
            }

            return;
        }
        else
        {
            std::cout<<"ENVI Error: Only BSQ operations are currently supported."<<std::endl;
            return;
        }
    }

    void getWavelength(void* band, double wavelength)
    {
        //this function retrieves a band via interpolation of wavelength values

        //get the low and high band numbers
        unsigned int low, high;
        get_surrounding_bands(wavelength, low, high);
        //std::cout<<"High: "<<high<<"   Low: "<<low<<std::endl;

		//calculate the position between the two bands
		double a;
		if(low == high)
			a = 1.0;
		else
			a = (wavelength - header.wavelength[low]) / (header.wavelength[high] - header.wavelength[low]);

        //read both bands
        void* A;
        A = malloc( getBandSize() );
        getBand(A, low);

        void* B;
        B = malloc( getBandSize() );
        getBand(B, high);

        //interpolate between the bands
        interpolate_band(band, A, B, a);

        //free the surrounding bands
        free(A);
        free(B);
    }

	//close file
	void close()
	{
        //std::cout<<"ENVI File Mode: "<<mode<<std::endl;

		//close the data stream
		fclose(data);

		//close the header
		if(mode != "r")
			header.save();
	}

};		//end EnviFile


#endif

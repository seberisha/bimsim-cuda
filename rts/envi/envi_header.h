#ifndef ENVI_HEADER_H
#define ENVI_HEADER_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>

//information from an ENVI header file
//A good resource can be found here: http://www.exelisvis.com/docs/enviheaderfiles.html

namespace Envi
{
    enum dataType {dummy0, int8, int16, int32, float32, float64, complex32, dummy7, dummy8, complex64, dummy10, dummy11, uint16, uint32, int64, uint64};
    enum interleaveType {BIP, BIL, BSQ};	//bip = Z,X,Y; bil = X,Z,Y; bsq = X,Y,Z
    enum endianType {endianLittle, endianBig};
}
struct EnviHeader
{
	std::string name;

	std::string description;

	unsigned int samples;	//x-axis
	unsigned int lines;		//y-axis
	unsigned int bands;		//spectral axis
	unsigned int header_offset;		//header offset for binary file (in bytes)
	std::string file_type;			//should be "ENVI Standard"

	Envi::dataType data_type;			//data representation; common value is 4 (32-bit float)


	Envi::interleaveType interleave;

	std::string sensor_type;		//not really used

	Envi::endianType byte_order;			//little = least significant bit first (most systems)

	double x_start, y_start;		//coordinates of the upper-left corner of the image
	std::string wavelength_units;	//stored wavelength units
	std::string z_plot_titles[2];

	double pixel_size[2];			//pixel size along X and Y

	std::vector<std::string> band_names;	//name for each band in the image
	std::vector<double> wavelength;		//wavelength for each band

	EnviHeader()
	{
        name = "";

		//specify default values for a new or empty ENVI file
		samples = 0;
		lines = 0;
		bands = 0;
		header_offset = 0;
		data_type = Envi::float32;
		interleave = Envi::BSQ;
		byte_order = Envi::endianLittle;
		x_start = y_start = 0;
		pixel_size[0] = pixel_size[1] = 1;

		//strings
		file_type = "ENVI Standard";
		sensor_type = "Unknown";
		wavelength_units = "Unknown";
		z_plot_titles[0] = z_plot_titles[1] = "Unknown";
	}

	std::string trim(std::string line)
	{
		//trims whitespace from the beginning and end of line
		unsigned int start_i, end_i;
		for(start_i=0; start_i < line.length(); start_i++)
			if(line[start_i] != 32)
			{
				break;
			}

		for(end_i = line.length()-1; end_i >= start_i; end_i--)
			if(line[end_i] != ' ')
			{
				break;
			}

		return line.substr(start_i, end_i - start_i+1);
	}


	std::string get_token(std::string line)
	{
		//returns a variable name; in this case we look for the '=' sign
		size_t i = line.find_first_of('=');

		std::string result;
		if(i != std::string::npos)
			result = trim(line.substr(0, i-1));

		return result;
	}

	std::string get_data_str(std::string line)
	{
		size_t i = line.find_first_of('=');

		std::string result;
		if(i != std::string::npos)
			result = trim(line.substr(i+1));
		else
		{
			std::cout<<"ENVI Header error - data not found for token: "<<get_token(line)<<std::endl;
			exit(1);
		}
		return result;
	}

	std::string get_brace_str(std::string token, std::string line, std::ifstream &file)
	{
		//this function assembles all of the characters between curly braces
		//this is how strings are defined within an ENVI file

		std::string result;

		//first, find the starting brace
		size_t i;
		do
		{
			i = line.find_first_of('{');
			if(i != std::string::npos)
				break;
		}while(file);

		//if i is still npos, we have reached the end of the file without a brace...something is wrong
		if(i == std::string::npos)
		{
			std::cout<<"ENVI Header error - string token declared without being defined: "<<token<<std::endl;
			exit(1);
		}
		line = line.substr(i+1);

		//copy character data into the result string until we find a closing brace
		while(file)
		{
			i = line.find_first_of('}');


			if(i != std::string::npos)
			{
				result += line.substr(0, i);
				break;
			}
			else
				result += line;

			getline(file, line);
		}

		if(i == std::string::npos)
		{
			std::cout<<"ENVI Header error - string token declared without a terminating '}': "<<token<<std::endl;
			exit(1);
		}

		return trim(result);
	}

	std::vector<std::string> get_string_seq(std::string token, std::string sequence)
	{
		//this function returns a sequence of comma-delimited strings
		std::vector<std::string> result;

		std::string entry;
		size_t i;
		do
		{
			i = sequence.find_first_of(',');
			entry = sequence.substr(0, i);
			sequence = sequence.substr(i+1);
			result.push_back(trim(entry));
		}while(i != std::string::npos);

		return result;
	}

	std::vector<double> get_double_seq(std::string token, std::string sequence)
	{
		//this function returns a sequence of comma-delimited strings
		std::vector<double> result;
		std::string entry;
		size_t i;
		do
		{
			i = sequence.find_first_of(',');
			entry = sequence.substr(0, i);
			sequence = sequence.substr(i+1);
			result.push_back(atof(entry.c_str()));
			//std::cout<<entry<<"   ";
		}while(i != std::string::npos);

		return result;
	}

	bool load(std::string filename)
	{
		//open the header file
		std::ifstream file(filename.c_str());

		if(!file)
		{
            return false;
		}

		//the first line should just be "ENVI"
		std::string line;
		file>>line;
		if(line != "ENVI")
		{
			std::cout<<"ENVI Header Error: The header doesn't appear to be an ENVI file. The first line should be 'ENVI'."<<std::endl;
			exit(1);
		}

		//for each line in the file, get the token
		std::string token;

        //get a line
        getline(file, line);
		while(file)
		{



			//get the token
			token = get_token(line);

			if(token == "description")
				description = get_brace_str(token, line, file);
			else if(token == "band names")
			{
				std::string string_sequence = get_brace_str(token, line, file);
				band_names = get_string_seq(token, string_sequence);
			}
			else if(token == "wavelength")
			{
				std::string string_sequence = get_brace_str(token, line, file);
				wavelength = get_double_seq(token, string_sequence);
			}
			else if(token == "pixel size")
			{
				std::string string_sequence = get_brace_str(token, line, file);
				std::vector<double> pxsize = get_double_seq(token, string_sequence);
				pixel_size[0] = pxsize[0];
				pixel_size[1] = pxsize[1];
			}
			else if(token == "z plot titles")
			{
				std::string string_sequence = get_brace_str(token, line, file);
				std::vector<std::string> titles = get_string_seq(token, string_sequence);
				z_plot_titles[0] = titles[0];
				z_plot_titles[1] = titles[1];
			}

			else if(token == "samples")
				samples = atoi(get_data_str(line).c_str());
			else if(token == "lines")
				lines = atoi(get_data_str(line).c_str());
			else if(token == "bands")
				bands = atoi(get_data_str(line).c_str());
			else if(token == "header offset")
				header_offset = atoi(get_data_str(line).c_str());
			else if(token == "file type")
				file_type = get_data_str(line);
			else if(token == "data type")
				data_type = (Envi::dataType)atoi(get_data_str(line).c_str());
			else if(token == "interleave")
			{
				std::string interleave_str = get_data_str(line);
				if(interleave_str == "bip")
					interleave = Envi::BIP;
				else if(interleave_str == "bil")
					interleave = Envi::BIL;
				else if(interleave_str == "bsq")
					interleave = Envi::BSQ;
			}
			else if(token == "sensor type")
				sensor_type = get_data_str(line);
			else if(token == "byte order")
				byte_order = (Envi::endianType)atoi(get_data_str(line).c_str());
			else if(token == "x start")
				x_start = atof(get_data_str(line).c_str());
			else if(token == "y start")
				y_start = atof(get_data_str(line).c_str());
			else if(token == "wavelength units")
				wavelength_units = get_data_str(line);

            //get the next line
            getline(file, line);
		}

		//make sure the number of bands matches the number of wavelengths
		unsigned int wavelengths = wavelength.size();
		if(bands != wavelengths)
		{
            std::cout<<"ENVI Header Error -- Number of wavelengths doesn't match the number of bands.  Bands = "<<bands<<", Wavelengths = "<<wavelength.size()<<std::endl;
            exit(1);
        }

		//close the file
		file.close();

		//set the file name
		name = filename;

		return true;
	}

	void save(std::string filename)
	{
		//open a file
		std::ofstream outfile(filename.c_str());

		//write the ENVI type identifier
		outfile<<"ENVI"<<std::endl;

		//output all of the data
		outfile<<"description = {"<<std::endl;
		outfile<<"  "<<description<<"}"<<std::endl;

		outfile<<"samples = "<<samples<<std::endl;
		outfile<<"lines = "<<lines<<std::endl;
		outfile<<"bands = "<<bands<<std::endl;
		outfile<<"header offset = "<<header_offset<<std::endl;
		outfile<<"file type = "<<file_type<<std::endl;
		outfile<<"data type = "<<data_type<<std::endl;
		outfile<<"interleave = ";
		if(interleave == Envi::BIP)
			outfile<<"bip";
		if(interleave == Envi::BIL)
			outfile<<"bil";
		if(interleave == Envi::BSQ)
			outfile<<"bsq";
			outfile<<std::endl;
		outfile<<"sensor type = "<<sensor_type<<std::endl;
		outfile<<"byte order = "<<byte_order<<std::endl;
		outfile<<"x start = "<<x_start<<std::endl;
		outfile<<"y start = "<<y_start<<std::endl;
		outfile<<"wavelength units = "<<wavelength_units<<std::endl;
		outfile<<"z plot titles = {";
			outfile<<z_plot_titles[0]<<", "<<z_plot_titles[1]<<"}"<<std::endl;
		outfile<<"pixel size = {"<<pixel_size[0]<<", "<<pixel_size[1]<<", units=Meters}"<<std::endl;
		if(band_names.size() > 0)
		{
			outfile<<"band names = {"<<std::endl;
			for(unsigned int i=0; i<band_names.size(); i++)
			{
				outfile<<band_names[i];
				if(i < band_names.size() - 1)
					outfile<<", ";
			}
			outfile<<"}"<<std::endl;
		}
		outfile<<"wavelength = {"<<std::endl;
			for(unsigned int i=0; i<wavelength.size()-1; i++)
				outfile<<wavelength[i]<<", ";
			outfile<<wavelength.back()<<"}"<<std::endl;

		outfile.close();
	}

	void save()
	{
        //std::cout<<"ENVI Header Name: "<<name<<std::endl;
		save(name);
	}

	//returns the size of a single value (in bytes)
	unsigned int valsize()
	{
		switch(data_type)
		{
		case Envi::int8:			//1 = 8-bit byte
			return 1;
		case Envi::int16:			//16-bit signed integer
		case Envi::uint16:		//16-bit unsigned integer
			return 2;
		case Envi::int32:			//32-bit signed long integer
		case Envi::float32:		//32-bit floating point
		case Envi::uint32:		//32-bit unsigned long integer
			return 4;
		case Envi::float64:		//64-bit double-precision floating point
		case Envi::complex32:		//32-bit complex value
		case Envi::int64:		    //64-bit signed long integer
		case Envi::uint64:		//64-bit unsigned long integer
			return 8;
		case Envi::complex64:		//64-bit complex value
			return 16;
		default:
			return 0;
		}

	}
};		//end EnviHeader

#endif


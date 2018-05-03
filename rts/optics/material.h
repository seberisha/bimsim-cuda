#ifndef MATERIALSTRUCT_H
#define MATERIALSTRUCT_H

#include <vector>
#include <ostream>
#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <sstream>
#include "rts/math/complex.h"
#include "rts/math/function.h"

//#define PI  3.14159

namespace rts{

enum field_type {field_microns, field_wavenumber, field_n, field_k, field_A, field_ignore};

//conversion functions

//convert wavenumber to lambda
template <class T>
static T _wn(T inverse_cm)
{
    return (T)10000.0/inverse_cm;
}

template <class T>
static T _2wn(T lambda)
{
    return (T)10000.0/lambda;
}

//convert absorbance to k
template <class T>
static T _A(T absorbance, T lambda)
{
    return (absorbance * lambda) / (4 * PI);
}
template <class T>
static T _2A(T k, T lambda)
{
    return (4 * PI * k)/lambda;
}

//define the dispersion as a single wavelength/refractive index pair
template <class T>
struct refIndex
{
    //wavelength (in microns)
    T lambda;
    complex<T> n;
};

template <class T>
struct entryType
{
    //list of value types per entry
    std::vector<field_type> valueList;

    entryType(std::string format)
    {
        //location of the end of a parameter
        size_t e;

        //string storing a token
        std::string token;

        do
        {
            //find the end of the first parameter
            e = format.find_first_of(',');

            //get the substring up to the comma
            token = format.substr(0, e);

            //turn the token into a val_type
            if(token == "microns")
                valueList.push_back(field_microns);
            else if(token == "wavenumber")
                valueList.push_back(field_wavenumber);
            else if(token == "n")
                valueList.push_back(field_n);
            else if(token == "k")
                valueList.push_back(field_k);
            else if(token == "A")
                valueList.push_back(field_A);
            else
                valueList.push_back(field_ignore);

            //remove the first token from the format string
            format = format.substr(e+1, format.length()-1);
        }while(e != std::string::npos);



    }

    void addValue(field_type value)
    {
        valueList.push_back(value);
    }

    refIndex<T> inputEntry(std::string line, T scaleA = 1.0)
    {
        T val;
        std::stringstream ss(line);

        //create a new refractive index
        refIndex<T> newRI;


        //read the entry from an input string
        for(unsigned int i=0; i<valueList.size(); i++)
        {


            while(ss.peek() < '0' || ss.peek() > '9')
            {
                //cout<<"bad char"<<endl;
                ss.ignore();
            }

            //retrieve the value
            ss>>val;
            //cout<<val<<endl;
            // std::cout<<"\t\t val: "<<val<<endl;

            //store the value in the appropriate location
            switch(valueList[i])
            {
            case field_microns:
                newRI.lambda = val;
                break;
            case field_wavenumber:
                newRI.lambda = _wn(val);
                break;
            case field_n:
                newRI.n.real(val);
                break;
            case field_k:
                newRI.n.imag(val);
                break;
            case field_A:
                newRI.n.imag(_A(val * scaleA, newRI.lambda));
                break;
            }
        }

        //return the refractive index associated with the entry
        return newRI;

    }

    std::string outputEntry(refIndex<T> material)
    {
        //std::string result;
        std::stringstream ss;

        //for each field in the entry
        for(int i=0; i<valueList.size(); i++)
        {
            if(i > 0)
                ss<<"\t";
            //store the value in the appropriate location
            switch(valueList[i])
            {
            case field_microns:
                ss<<material.lambda;
                break;
            case field_wavenumber:
                ss<<_2wn(material.lambda);
                break;
            case field_n:
                ss<<material.n.real();
                break;
            case field_k:
                ss<<material.n.imag();
                break;
            case field_A:
                ss<<_2A(material.n.imag(), material.lambda);
                break;
            }

        }
        return ss.str();
    }
};


//a material is a list of refractive index values
template <class T>
class material
{
    std::string name;
    //dispersion (refractive index as a function of wavelength)
    //std::vector< refIndex<T> > dispersion;
    function< T, complex<T> > dispersion;

    //average refractive index (approximately 1.4)
    T n0;

    /*void add(refIndex<T> ri)
        {
            //refIndex<T> converted = convert(ri, measurement);
            dispersion.push_back(ri);
        }*/

    //comparison function for sorting
    static bool compare(refIndex<T> a, refIndex<T> b)
    {
        return (a.lambda < b.lambda);
    }

    //comparison function for searching lambda
    /*static bool findCeiling(refIndex<T> a, refIndex<T> b)
        {
            return (a.lambda > b.lambda);
        }*/

public:
    void add(T lambda, complex<T> n)
    {
      //  std::cout<<"\t\t====>in add "<<lambda<<","<<n<<std::endl;
        dispersion.insert(lambda, n);
    }

    std::string getName()
    {
        return name;
    }
    void setName(std::string n)
    {
        name = n;
    }
    T getN0()
    {
        return n0;
    }
    void setN0(T n)
    {
        n0 = n;
    }

    void setM(function< T, complex<T> > m)
    {
        dispersion = m;
    }
    unsigned int nSamples()
    {
        //return dispersion.size();
        //sb: size wasn't part of function.h in rts...guessing getN is the same
        return dispersion.getN();
    }
    material<T> computeN(T _n0, unsigned int n_samples = 0, T pf = 2)
    {
        /*	This function computes the real part of the refractive index
                                from the imaginary part.  The Hilbert transform is required. I
                                use an FFT in order to simplify this, so either the FFTW or CUFFT
                                packages are required.  CUFFT is used if this file is passed through
                                a CUDA compiler.  Otherwise, FFTW is used if available.
                        */
        n0 = _n0;

        int N;
        if(n_samples)
            N = n_samples;
        else
            //N = dispersion.size();
            //sb: no size() in function.h of rts...replacing it with getN()
            N = dispersion.getN();

#ifdef FFTW_AVAILABLE
        //allocate memory for the FFT
        complex<T>* Chi2 = (complex<T>*)fftw_malloc(sizeof(complex<T>) * N * pf);
        complex<T>* Chi2FFT = (complex<T>*)fftw_malloc(sizeof(complex<T>) * N * pf);
        complex<T>* Chi1 = (complex<T>*)fftw_malloc(sizeof(complex<T>) * N * pf);

        //create an FFT plan for the forward and inverse transforms
        fftw_plan planForward, planInverse;
        planForward = fftw_plan_dft_1d(N*pf, (fftw_complex*)Chi2, (fftw_complex*)Chi2FFT, FFTW_FORWARD, FFTW_ESTIMATE);
        planInverse = fftw_plan_dft_1d(N*pf, (fftw_complex*)Chi2FFT, (fftw_complex*)Chi1, FFTW_BACKWARD, FFTW_ESTIMATE);

        float k, alpha;
        T chi_temp;

        //the spectrum will be re-sampled in uniform values of wavenumber
//        T nuMin = _2wn(dispersion.back().lambda);
//        T nuMax = _2wn(dispersion.front().lambda);


        //sb: added these to replace the above 2 lines
        T nuMin = _2wn(dispersion.backLambda());
        T nuMax = _2wn(dispersion.frontLambda());

        std::cout<<"\t\t nuMin ---"<<nuMin<<std::endl;
        std::cout<<"\t\t nuMax ---"<<nuMax<<std::endl;

        T dnu = (nuMax - nuMin)/(N-1);
        T lambda, tlambda;
        for(int i=0; i<N; i++)
        {
            //go from back-to-front (wavenumber is the inverse of wavelength)
            lambda = _wn(nuMax - i * dnu);

            //compute the frequency
            k = 2 * PI / (lambda);

            //get the absorbance
            alpha = getN(lambda).imag() * (2 * k);

            //compute chi2
            Chi2[i] = -alpha * (n0 / k);
        }



        //use linear interpolation between the start and end points to pad the spectrum
        //complex<T> nMin = dispersion.back();
        //complex<T> nMax = dispersion.front();
        T a;
        for(int i=N; i<N*pf; i++)
        {
            //a = (T)(i-N)/(T)(N*pf - N);
            //Chi2[i] = a * Chi2[0] + ((T)1 - a) * Chi2[N-1];

            Chi2[i] = 0.0;//Chi2[N-1];
        }



        //perform the FFT
        fftw_execute(planForward);


        //perform the Hilbert transform in the Fourier domain
        complex<T> j(0, 1);
        for(int i=0; i<N*pf; i++)
        {
            //if w = 0, set the DC component to zero
            if(i == 0)
                Chi2FFT[i] *= (T)0.0;
            //if w <0, multiply by i
            else if(i < N*pf/2.0)
                Chi2FFT[i] *= j;
            //if i > N/2, multiply by -i
            else
                Chi2FFT[i] *= -j;
        }

        //execute the inverse Fourier transform (completing the Hilbert transform)
        fftw_execute(planInverse);



        //divide the Chi1 values by N
        for(int i=0; i<N*pf; i++)
            Chi1[i] /= (T)(N*pf);

        //create a new material
        material<T> newM;

        newM.dispersion.clear();
       // std::cout<<"\t\t after newM.dispersion.clear==="<<std::endl;
        refIndex<T> ri;
        for(int i=0; i<N; i++)
        {
            ri.lambda = _wn(nuMax - i * dnu);
            ri.n.real(Chi1[i].real() / (2 * n0) + n0);
            ri.n.imag(getN(ri.lambda).imag());

            //newM.dispersion.push_back(ri);
            //sb: use insert instead

           // std::cout<<"\t\t ri.lambda==="<<ri.lambda<<std::endl;
//            std::cout<<"\t\t ri.n.real()==="<<ri.n.real()<<std::endl;
//            std::cout<<"\t\t ri.n.imag()==="<<ri.n.imag()<<std::endl;

            newM.dispersion.insert(ri.lambda, ri.n);
        }





        //dispersion[i].n.real(Chi1[i].real() / (2 * n0) + n0);


        /*//output the Absorbance value
                        ofstream outOrig("origN.txt");
                        for(int i=0; i<N; i++)
                                outOrig<<dispersion[i].lambda<<"     "<<dispersion[i].n.real()<<endl;

                        //output the Chi2 value
                        ofstream outChi2("chi2.txt");
                        for(int i=0; i<N; i++)
                                outChi2<<dispersion[i].lambda<<"     "<<Chi2[i].real()<<endl;

                        //output the Fourier transform
                        ofstream outFFT("chi2_FFT.txt");
                        for(int i=0; i<N; i++)
                        {
                        float mag = std::sqrt( std::pow(Chi2FFT[i].real(), 2.0) + std::pow(Chi2FFT[i].imag(), 2.0));
                                outFFT<<dispersion[i].lambda<<"     "<<mag<<endl;
                        }

                        //output the computed Chi1 value
                        ofstream outChi1("chi1.txt");
                        for(int i=0; i<N; i++)
                        {
                                outChi1<<dispersion[i].lambda<<"     "<<Chi1[i].real()<<"     "<<Chi1[i].imag()<<endl;
                        }

                        ofstream outN("n.txt");
                        for(int i=0; i<N; i++)
                                outN<<dispersion[i].lambda<<"     "<<Chi1[i].real() / (2 * n0) + n0<<endl;*/


        //de-allocate memory
        fftw_destroy_plan(planForward);
        fftw_destroy_plan(planInverse);
        fftw_free(Chi2);
        fftw_free(Chi2FFT);
        fftw_free(Chi1);

        return newM;
#else
        std::cout<<"MATERIAL Error: Must have FFTW in order to compute Kramers-Kronig."<<std::endl;
        return material<T>();
#endif
    }

    material(T lambda = 1.0, T n = 1.4, T k = 0.0)
    {
        dispersion.insert(lambda, complex<T>(0.0, k));
        /*//create a default refractive index
            refIndex<T> def;
            def.lambda = lambda;
            def.n.real(n);
            def.n.imag(k);
            add(def);
                        */
        //set n0
        n0 = n;
    }

    material(std::string filename, std::string format = "", T scaleA = 1.0)
    {
        //std::cout<<"\t\t\t material"<<std::endl;

        fromFile(filename, format);
        n0 = 0;
    }

    void fromFile(std::string filename, std::string format = "", T scaleA = 1.0)
    {
        name = filename;
        //clear any previous values
        dispersion = rts::function< T, complex<T> >();

        //printf("\t\t filename: %s\n",filename.c_str());

        //load the file into a string
        std::ifstream ifs(filename.c_str());

        std::string line;

        if(!ifs.is_open())
        {
            std::cout<<"Material Error -- file not found: "<<filename<<std::endl;
            exit(1);
        }

        //process the file as a string
        std::string instr((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
        fromStr(instr, format, scaleA);
    }

    void fromStr(std::string str, std::string format = "", T scaleA = 1.0)
    {

        //std::cout<<"\t\t\t============>fromStr"<<std::endl;

        //create a string stream to process the input data
        std::stringstream ss(str);

        //this string will be read line-by-line (where each line is an entry)
        std::string line;

        //if the format is not provided, see if it is in the file, otherwise use a default
        if(format == "")
        {
            //set a default format of "lambda,n,k"
            format = "microns,n,k";

            //see if the first line is a comment
            char c = ss.peek();
            if(c == '#')
            {
              //  std::cout<<"\t\t\t============>if(c == '#')"<<std::endl;

                //get the first line
                getline(ss, line);
                //look for a bracket, denoting the start of a format string
                int istart = line.find('[');
                if(istart != std::string::npos)
                {
                    //look for a bracket terminating the format string
                    int iend = line.find(']');
                    if(iend != std::string::npos)
                    {
                        //read the string between the brackets
                        format = line.substr(istart+1, iend - istart - 1);
                    }
                }
            }

        }

        entryType<T> entry(format);

     //   std::cout<<"Loading material with format: "<<format<<std::endl;

        while(!ss.eof())
        {
            //read a line from the string
            getline(ss, line);
          //  printf("\t\t line: %s\n",line.c_str());
            //if the line is not a comment, process it
            if(line[0] != '#')
            {
                //load the entry and add it to the dispersion list
                add(entry.inputEntry(line, scaleA).lambda, entry.inputEntry(line, scaleA).n);
            }
            //generally have to peek to trigger the eof flag
            ss.peek();
        }

        //sort the vector by lambda
        //sort(dispersion.begin(), dispersion.end(), &material<T>::compare);
    }

    //convert the material to a string
    std::string toStr(std::string format = "microns,n,k", bool reverse_order = false)
    {
        std::stringstream ss;
        entryType<T> entry(format);
        //create a new refractive index
        refIndex<T> ri;


        if(reverse_order)
        {
//            for(int l=dispersion.size() - 1; l>=0; l--)
//            {
//                if(l < dispersion.size() - 1) ss<<std::endl;
//                ss<<entry.outputEntry(dispersion[l]);
//            }

            //sb: no size() in function.h of rts...replacing it with getN()
            for(int l=dispersion.getN() - 1; l>=0; l--)
            {

                ri.lambda = dispersion[l].x;
                ri.n = dispersion[l].y;
                if(l < dispersion.getN() - 1) ss<<std::endl;
                //ss<<entry.outputEntry(dispersion[l]);
                ss<<entry.outputEntry(ri);
            }

        }
        else
        {
            /*
            for(unsigned int l=0; l<dispersion.size(); l++)
            {
                if(l > 0) ss<<std::endl;
                ss<<entry.outputEntry(dispersion[l]);
            }*/

            //sb: no size() in function.h of rts...replacing it with getN()


            for(unsigned int l=0; l<dispersion.getN(); l++)
            {
                ri.lambda = dispersion[l].x;
                ri.n = dispersion[l].y;
                if(l > 0) ss<<std::endl;
                //ss<<entry.outputEntry(dispersion[l]);
                ss<<entry.outputEntry(ri);
            }
        }

        return ss.str();
    }

    void save(std::string filename, std::string format = "microns,n,k", bool reverse_order = false)
    {
        std::ofstream outfile(filename.c_str());
        outfile<<"#material file saved as [" + format + "]"<<std::endl;
        outfile<<toStr(format, reverse_order)<<std::endl;

    }

    //convert between wavelength and wavenumber
    /*void nu2lambda(T s = (T)1)
        {
            for(int i=0; i<dispersion.size(); i++)
                dispersion[i].lambda = s/dispersion[i].lambda;
        }

        void lambda2nu(T s = (T)1)
        {
            for(int i=0; i<dispersion.size(); i++)
                dispersion[i].lambda = s/dispersion[i].lambda;
        }*/


    refIndex<T>& operator[](unsigned int i)
    {
        return dispersion[i];

    }

    complex<T> getN(T l)
    {
        //return complex<T>(1.0, 0.0);
        complex<T> ri = dispersion.linear(l) + n0;
        return ri;
    }

    function<T, complex<T> > getF()
    {
        return dispersion + complex<T>(n0, 0.0);
    }

    //returns the scattering efficiency at wavelength l
    complex<T> eta(T l)
    {
        //get the refractive index
        complex<T> ri = getN(l);

        //convert the refractive index to scattering efficiency
        return ri*ri - 1.0;

    }
    //interpolate the given lambda value and return the index of refraction
    complex<T> operator()(T l)
    {
        //std::cout<<"\t\t\t operator()(T l)"<<std::endl;
        return getN(l);
    }


};
}   //end namespace rts

template <typename T>
std::ostream& operator<<(std::ostream& os, rts::material<T> m)
{
    os<<m.toStr();

    return os;
}



#endif

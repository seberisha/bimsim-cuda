#ifndef RTS_FUNCTION_H
#define RTS_FUNCTION_H

#include <string>

namespace rts{

template <class X, class Y>
class function
{
    //datapoint class for storing function points
    struct dataPoint
    {
        X x;
        Y y;
    };

    //function data
    std::vector<dataPoint> f;

    //comparison function for searching lambda
    static bool findCeiling(dataPoint a, dataPoint b)
    {
        return (a.x > b.x);
    }


public:

    //adding these functions that were missing from previous rts:function.h code
    void clear()
    {
        f.clear();
    }

    X backLambda()
    {
        return f.back().x;
    }

    X frontLambda()
    {
        return f.front().x;
    }

    Y linear(X x)
    {
        //declare an iterator
        typename std::vector< dataPoint >::iterator it;

        dataPoint s;
        s.x = x;

        it = search(f.begin(), f.end(), &s, &s + 1, &function<X, Y>::findCeiling);

        //if the wavelength is past the end of the list, return the back
        if(it == f.end())
            return f.back().y;
        //if the wavelength is before the beginning of the list, return the front
        else if(it == f.begin())
            return f.front().y;
        //otherwise interpolate
        else
        {
            X xMax = (*it).x;
            X xMin = (*(it - 1)).x;
            //std::cout<<lMin<<"----------"<<lMax<<std::endl;

            X a = (x - xMin) / (xMax - xMin);
            Y riMin = (*(it - 1)).y;
            Y riMax = (*it).y;
            Y interp;
            interp = riMin * a + riMax * (1.0 - a);
            return interp;
        }
    }

    void insert(X x, Y y)
    {
        //declare an iterator
        typename std::vector< dataPoint >::iterator it;

        dataPoint s;
        s.x = x;
        s.y = y;

        it = search(f.begin(), f.end(), &s, &s + 1, &function<X, Y>::findCeiling);

        //if the function value is past the end of the vector, add it to the back
        if(it == f.end())
            return f.push_back(s);
        //otherwise add the value at the iterator position
        else
        {
            f.insert(it, s);
        }

    }

    X getX(unsigned int i)
    {
        return f[i].x;
    }

    Y getY(unsigned int i)
    {
        return f[i].y;
    }

    unsigned int getN()
    {
        return f.size();
    }

    dataPoint operator[](int i)
    {
        return f[i];
    }

    function<X, Y> operator+(Y r)
    {
        function<X, Y> result;

        //add r to every point in f
        for(int i=0; i<f.size(); i++)
        {
            result.f.push_back(f[i]);
            result.f[i].y += r;
        }

        return result;
    }


};

}	//end namespace rts


#endif

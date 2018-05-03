#ifndef __PLANEWAVE__
#define __PLANEWAVE__

#include <iostream>
#include <sstream>

#include "rts/math/vector.h"

template<class T>
class planewave
{
    rts::vector<T, 3> k;
    rts::vector<T, 3> E;

    public:

    //constructor, initialize to an x-polarized wave propagating along z
    planewave()
    {
        k = rts::vector<T, 3>(0, 0, 1);
        E = rts::vector<T, 3>(1, 0, 0);
    }

    planewave(rts::vector<T, 3> k_vec, rts::vector<T, 3> E0)
    {
        k = k_vec;

        //enforce k \dot E = 0
        rts::vector<T, 3> s = E0.cross(k);
        rts::vector<T, 3> E_hat;

        if(s.len() == 0)
            E_hat = rts::vector<T, 3>(0, 0, 0);
        else
            E_hat = (s.cross(k)).norm();

        E = E_hat * (E_hat.dot(E0));
    }

    // refract will bend the wave vector k to correspond to the normalized vector v
    void refract(rts::vector<T, 3> v)
    {
        //make sure that v is normalized
        v = v.norm();


    }

    std::string toStr()
	{
		std::stringstream ss;

		ss<<"k = "<<k<<std::endl;
		ss<<"E = "<<E<<std::endl;

		return ss.str();
	}



};

template <typename T>
std::ostream& operator<<(std::ostream& os, planewave<T> p)
{
    os<<p.toStr();
    return os;
}



#endif

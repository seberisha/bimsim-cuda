#ifndef RTS_MATRIX_H
#define RTS_MATRIX_H

//#include "rts/vector.h"
#include <string.h>
#include <iostream>

namespace rts
{

template <class T, int N>
struct matrix
{
	//the matrix will be stored in column-major order (compatible with OpenGL)
	T M[N*N];

	matrix()
	{
		for(int r=0; r<N; r++)
			for(int c=0; c<N; c++)
				if(r == c)
					(*this)(r, c) = 1;
				else
					(*this)(r, c) = 0;
	}

	T& operator()(int row, int col)
	{
		return M[col * N + row];
	}

	matrix<T, N> operator=(T rhs)
	{
		int Nsq = N*N;
		for(int i=0; i<Nsq; i++)
			M[i] = rhs;

		return *this;
	}

	/*matrix<T, N> operator=(matrix<T, N> rhs)
	{
		for(int i=0; i<N; i++)
			M[i] = rhs.M[i];

		return *this;
	}*/

	vector<T, N> operator*(vector<T, N> rhs)
	{
		vector<T, N> result;

		for(int r=0; r<N; r++)
			for(int c=0; c<N; c++)
				result[r] += (*this)(r, c) * rhs[c];

		return result;
	}

	std::string toStr()
	{
		std::stringstream ss;

		for(int r = 0; r < N; r++)
		{
			ss<<"| ";
			for(int c=0; c<N; c++)
			{
				ss<<(*this)(r, c)<<" ";
			}
			ss<<"|"<<std::endl;
		}

		return ss.str();
	}




};

}	//end namespace rts

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, rts::matrix<T, N> M)
{
    os<<M.toStr();
    return os;
}

//#if __GNUC__ > 3 && __GNUC_MINOR__ > 7
//template<class T, int N> using rtsMatrix = rts::matrix<T, N>;
//#endif

#endif

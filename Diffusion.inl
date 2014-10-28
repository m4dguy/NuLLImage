#ifndef DIFFUSION_INL
#define DIFFUSION_INL

#include "Diffusion.h"

//linear diffusion
//blurs the image; related to Gauss filtering
template <typename T> void Diffusion::linearDiffusion(const Matrix<T>& src, Matrix<T>& dst, const double stepsize, const uint iterations)
{
	T val;
	const size_t width = src.width();
	const size_t height = src.height();
	
	Matrix<T> tmp(width, height);
	dst = src;

	//main loop; steps derived from jacobi-iteration
	for(uint i=0; i<iterations; ++i)
	{
		for(uint y=0; y<height; ++y)
		{
			for(uint x=0; x<width; ++x)
			{
				val = dst(x,y) + stepsize * (dst.getMirrored(x+1,y) - 2*(dst(x,y)) + dst.getMirrored(x-1,y) + dst.getMirrored(x,y+1) - (2*dst(x,y)) + dst.getMirrored(x,y-1));
				tmp(x,y) = val;
			}
		}
		dst.swap(tmp);
	}
}

#endif
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

	//main loop: steps derived with forward time difference scheme
	for(uint i=0; i<iterations; ++i)
	{
		for(uint y=0; y<height; ++y)
		{
			for(uint x=0; x<width; ++x)
			{
				val = dst.getMirrored(x + 1, y) + dst.getMirrored(x - 1, y) + dst.getMirrored(x, y + 1) + dst.getMirrored(x, y - 1) - 4 * (dst(x, y));
				tmp(x,y) = dst(x,y) + stepsize * val;
			}
		}
		dst.swap(tmp);
	}
}

//nonlinear diffusion with with perona-malik diffusity
//implementation based on "Scale-Space and Edge Detection Using Anisotropic Diffusion" by P.Perona and J.Malik
template <typename T> void Diffusion::nonlinearDiffusion(const Matrix<T>& src, Matrix<T>& dst, const double lambda, const double stepsize, const uint iterations)
{
	const double lambdaSqInv = 1. / (lambda * lambda);
	const size_t width = src.width();
	const size_t height = src.height();

	Matrix<T> tmp(width, height);
	dst = src;

	T val;
	T cN, cS, cE, cW;		//edge stopping values
	T dN, dS, dE, dW;		//one-sided differences

	//main loop: steps derived with forward time difference scheme
	for (uint i = 0; i<iterations; ++i)
	{
		for (uint y = 0; y<height; ++y)
		{
			for (uint x = 0; x<width; ++x)
			{
				dN = dst.getMirrored(x, y - 1) - dst(x, y);
				dS = dst.getMirrored(x, y + 1) - dst(x, y);
				dE = dst.getMirrored(x - 1, y) - dst(x, y);
				dW = dst.getMirrored(x + 1, y) - dst(x, y);

				cN = 1. / (1. + (dN*dN * lambdaSqInv));
				cS = 1. / (1. + (dS*dS * lambdaSqInv));
				cE = 1. / (1. + (dE*dE * lambdaSqInv));
				cW = 1. / (1. + (dW*dW * lambdaSqInv));

				val = cN * dN + cS * dS + cE * dE + cW * dW;

				tmp(x, y) = dst(x,y) + stepsize * val;
			}
		}
		dst.swap(tmp);
	}
}

#endif
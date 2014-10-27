#ifndef VARIATION_INL
#define VARIATION_INL

#include "Variation.h"

//variational blur
//removes noise by blurring while trying to preserve edges
template <typename T> void Variation::variationalBlur(const Matrix<T>& src, Matrix<T>& dst, const double alpha, const uint iterations)
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
				val = src(x,y) + alpha * (dst.getMirrored(x-1,y) + dst.getMirrored(x+1,y) + dst.getMirrored(x,y-1) + dst.getMirrored(x,y+1));
				val /= (1 + alpha * 4);
				tmp(x,y) = val;
			}
		}
		dst.swap(tmp);
	}
}

//TODO
//variational blur using charbonnier diffusity/ edge stopping term
/*template <typename T> void Variation::charbonnierDenoise(const Matrix<T>& src, Matrix<T>& dst, const double alpha, const double lambda, const uint iterations)
{
	T val;
	const size_t width = src.width();
	const size_t height = src.height();
	
	Matrix<T> tmp(width, height);
	dst = src;

	const double tau = 0.1;

	const T lambdaSq = lambda * lambda;
	T beta, sqrtBeta, sqrtBeta3, grad;
	T diag, offdiag;
	T dx, dy, dxx, dxy, dyy;
	T a, b;

	//main loop; steps derived from jacobi-iteration
	for(uint i=0; i<iterations; ++i)
	{
		for(uint y=0; y<height; ++y)
		{
			for(uint x=0; x<width; ++x)
			{
				dx = (dst.getMirrored(x+1,y) - dst.getMirrored(x-1,y)) / 2.;
				dy = (dst.getMirrored(x,y+1) - dst.getMirrored(x,y-1)) / 2.;
				dxx = dst.getMirrored(x+1,y) - 2.*dst(x,y) + dst.getMirrored(x-1,y);
				dyy = dst.getMirrored(x,y+1) - 2.*dst(x,y) + dst.getMirrored(x,y-1);

				dxy = .25*((dst.getMirrored(x+1,y+1) - dst.getMirrored(x+1,y-1)) - (dst.getMirrored(x-1,y+1) - dst.getMirrored(x-1,y-1)));
				dxy = (dst.getMirrored(x-1,y) + dst.getMirrored(x+1,y) + dst.getMirrored(x,y-1) + dst.getMirrored(x,y+1) - 2*dst(x,y) - dst.getMirrored(x-1,y+1) + dst.getMirrored(x+1,y-1)) / 2.;

				beta = 1. + (dx*dx+dy*dy) / (lambdaSq);
				sqrtBeta = sqrt(beta);
				sqrtBeta3 = sqrtBeta * sqrtBeta * sqrtBeta;

				//diag = 1. - ((alpha * 4.) / sqrtBeta);

				//a = (dxx+dxy) / (2 * lambdaSq);
				//b = (dyy+dxy) / (2 * lambdaSq);

				//val = src(x,y) - (alpha * (dst.getMirrored(x+1,y)*(beta+a) + dst.getMirrored(x-1,y)*(beta-a) + dst.getMirrored(x,y+1)*(beta+b) + dst.getMirrored(x,y-1)*(beta-b))/sqrtBeta3);
				//tmp(x,y) = val / diag;


				grad = 4*dst(x,y)-(dst.getMirrored(x+1,y) + dst.getMirrored(x-1,y) + dst.getMirrored(x,y+1) + dst.getMirrored(x,y-1));
				
				val = alpha*((dx*(dxx+dxy)+dy*(dyy+dxy))/lambdaSq - (grad * beta));
				val /= sqrtBeta3;
				tmp(x,y) = dst(x,y) + tau * ((dst(x,y) - src(x,y)) - val);
			}
		}
		dst.swap(tmp);
	}
}*/


#endif
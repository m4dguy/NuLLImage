#ifndef NULLCONVOLVE_INL
#define NULLCONVOLVE_INL

#include "NuLLConvolve.h"

//discrete convolution of matrix (2D convolution)
template<typename T> void NuLLConvolve::convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Matrix<T>& kernel)
{
    T val;
    size_t width = mtx.width();
    size_t height = mtx.height();

    size_t kWidth = kernel.width();
    size_t kHeight = kernel.height();

    size_t offsetX = kernel.width() / 2;
    size_t offsetY = kernel.height() / 2;

    for(uint y=0; y<height; ++y)
    {
        for(uint x=0; x<width; ++x)
        {
            val = 0.0f;
            for(uint ky=0; ky<kHeight; ++ky)
            {
                for(uint kx=0; kx<kWidth; ++kx)
                {
                    val += kernel(kx, ky) * mtx.getMirrored(x+(offsetX-kx), y+(offsetY-ky));
                }
            }
            dst(x,y) = val;
        }
    }
}

//covolution with 1D kernel
//faster than 2D convolution, with linear running time
//valid alternative if kernel is separable
template<typename T> void NuLLConvolve::convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernel)
{
	convolve(mtx, dst, kernel, kernel);
}

//covolution with 1D kernel
//faster than 2D convolution, with linear running time
//valid alternative if kernel is separable
template <typename T> void NuLLConvolve::convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelX, const Vector<T>& kernelY)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
	convolveX(mtx, tmp, kernelX);
	convolveY(tmp, dst, kernelY);
}

//convolution with 1D kernel
//x direction only
template <typename T> void NuLLConvolve::convolveX(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelX)
{
	const size_t width = mtx.width();
	const size_t height = mtx.height();

	T val;
	const size_t kSizeX = kernelX.size();
	const size_t offsetX = kSizeX / 2;

	//convolution in x-direction
	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			val = 0.;
			for(uint i=0; i<kSizeX; ++i)
			{
				val += kernelX[i] * mtx.getMirrored(x+offsetX-i, y);
			}
			dst(x,y) = val;
		}
	}
}

//convolution with 1D kernel
//x direction only
template <typename T> void NuLLConvolve::convolveY(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelY)
{
	const size_t width = mtx.width();
	const size_t height = mtx.height();

	T val;
	const size_t kSizeY = kernelY.size();
	const size_t offsetY = kSizeY / 2;

	//convolution in y-direction
	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			val = 0;
			for(uint i=0; i<kSizeY; ++i)
			{
				val += kernelY[i] * mtx.getMirrored(x, y+offsetY-i);
			}
			dst(x,y) = val;
		}
	}
}

//normalizes matrix to make sum of entries equal 1
//useful for kernel normalization
template<typename T> void NuLLConvolve::normalize(Matrix<T>& mtx)
{
    T factor = 0;
    const size_t height = mtx.height();
    const size_t width = mtx.width();

    for(uint y=0; y<height; ++y)
        for(uint x=0; x<width; ++x)
            factor += mtx(x,y);

    mtx /= factor;
}

//creates identity kernel for convolution
//useful for kernel design
template <typename T> void NuLLConvolve::identityKernel(Matrix<T>& dst, const int radius)
{
	dst.fill(0.0);
	dst(radius, radius) = 1.0;
}

//disc-shaped pillbox kernel
template <typename T> void NuLLConvolve::pillboxKernel(Matrix<T>& dst, const int radius)
{
    const int dim = (2 * radius) + 1;
    int radSq = radius * radius;
    T scale = 0;

    //kernel calculation
    for(int y=0; y<dim; ++y)
    {
        for(int x=0; x<dim; ++x)
        {
            if(((radius-x)*(radius-x) + (radius-y) * (radius-y)) <= radSq)
            {
                dst(x,y) = 1.0;
                ++scale;
            }
        }
    }
    dst /= scale;
}

//creates a gaussian kernel of given size
//sigma/ variance parameter only has notable effect for large kernels
template <typename T> void NuLLConvolve::gaussianKernel(Matrix<T>& dst, const int radius, const double sigma)
{
	if(!sigma)
		sigma = radius;

    const uint dim = (2 * radius) + 1;
	const double pi = 2.0 * asinf(1.0);

    T div = 2 * pi * sigma * sigma;
    T kx, ky, res, scale;
    scale = 0;

    //kernel calculation
    for(uint y=0; y<dim; ++y)
    {
        for(uint x=0; x<dim; ++x)
        {
            kx = (x - radius) * (x - radius);
            ky = (y - radius) * (y - radius);
            res = -(ky + kx);
            res = exp(res / (2 * sigma * sigma));
            res /= div;
            dst(x,y) = res;

            scale += res;
        }
    }

	//normalization: sum of entries only roughly == 1
    dst /= scale;
}

#endif

#ifndef NULLPROCESSING_INL
#define NULLPROCESSING_INL

#include "NuLLProcessing.h"

//discrete convolution of matrix (2D convolution)
template<typename T> void NuLLProcessing::convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Matrix<T>& kernel)
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
template<typename T> void NuLLProcessing::convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernel)
{
	convolve(mtx, dst, kernel, kernel);
}

//covolution with 1D kernel
//faster than 2D convolution, with linear running time
//valid alternative if kernel is separable
template <typename T> void NuLLProcessing::convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelX, const Vector<T>& kernelY)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
	NuLLProcessing::convolveX(mtx, tmp, kernelX);
	NuLLProcessing::convolveY(tmp, dst, kernelY);
}

//convolution with 1D kernel
//x direction only
template <typename T> void NuLLProcessing::convolveX(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelX)
{
	const uint width = mtx.width();
	const uint height = mtx.height();

	T val;
	const uint kSizeX = kernelX.size();
	const uint offsetX = kSizeX / 2;

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
template <typename T> void NuLLProcessing::convolveY(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelY)
{
	const uint width = mtx.width();
	const uint height = mtx.height();

	T val;
	const uint kSizeY = kernelY.size();
	const uint offsetY = kSizeY / 2;

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
template<typename T> void NuLLProcessing::normalize(Matrix<T>& mtx)
{
    T factor = 0;
    size_t height = mtx.height();
    size_t width = mtx.width();

    for(uint y=0; y<height; ++y)
        for(uint x=0; x<width; ++x)
            factor += mtx(x,y);

    mtx /= factor;
}

//creates identity kernel for convolution
//useful for kernel design
template <typename T> void NuLLProcessing::identityKernel(Matrix<T>& dst, int radius)
{
	dst.fill(0.0);
	dst(radius, radius) = 1.0;
}

//disc-shaped pillbox kernel
template <typename T> void NuLLProcessing::pillboxKernel(Matrix<T>& dst, int radius)
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
template <typename T> void NuLLProcessing::gaussianKernel(Matrix<T>& dst, int radius, double sigma)
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

//smoothing with pillbox kernel
//uses separation theorem for fast convolution/ linear running time
template <typename T> void NuLLProcessing::boxBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	uint kSize = (2 * radius) + 1;
	Vector<T> kernel(kSize);
	kernel.fill(1.);
	kernel /= (T)kSize;

	convolve(mtx, dst, kernel);
}

//gaussian smoothing for matrices
//uses separation theorem for fast convolution/ linear running time
template <typename T> void NuLLProcessing::gaussianBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius, double sigma)
{
	if(!sigma)
		sigma = radius;

	const int kSize = (2 * radius) + 1;
	const double pi = 2.0 * asinf(1.0);
	const double divExp = 2 * sigma * sigma;
	const double div = sqrt(2 * pi) * sigma;
    T res, scale;
    scale = 0;

	Vector<T> kernel(kSize);
	for(int i=0; i<kSize; ++i)
	{
		res = exp(- ((i-radius)*(i-radius)) / divExp);
		res /= div;
		kernel[i] = res;
		scale += res;
	}
	kernel /= (T)scale;

	convolve(mtx, dst, kernel);
}

//apply affine rescale afterwards or else over and undershoots occur
template <typename T> void NuLLProcessing::differenceOfGaussians(const Matrix<T>& mtx, Matrix<T>& dst, int sigma1, int sigma2)
{
	Matrix<T> tmp(mtx.width(), mtx.height());

	if(sigma1)
		gaussianBlur(mtx, dst, sigma1);
	else
		dst = mtx;

	gaussianBlur(mtx, tmp, sigma2);
	dst -= tmp;
}

//first order derivative (central differential quotient) for matrices
template <typename T> void NuLLProcessing::firstDerivative(const Matrix<T>& mtx, Matrix<T>& dstGrad, Matrix<T>& dstDir)
{
    uint width = mtx.width();
    uint height = mtx.height();
    T gradMag;

	Matrix<T> dx(width, height);
	Matrix<T> dy(width, height);

	Vector<T> stencil(3);
	stencil[0] = -.5;
	stencil[1] = 0;
	stencil[2] = .5;

	convolveX(mtx, dx, stencil);
	convolveY(mtx, dy, stencil);

    for(uint y=0; y<height; ++y)
    {
        for(uint x=0; x<width; ++x)
        {
            gradMag = sqrt(dx(x,y) * dx(x,y) + dy(x,y) * dy(x,y));
            dstGrad(x,y) = gradMag;
			dstDir(x,y) = (T)atan2(dy(x,y), dx(x,y));
        }
    }
}

//hack: overloaded function for compatibility
//fix this...
template <typename T> void NuLLProcessing::firstDerivative(const Matrix<T>& mtx, Matrix<T>& dst)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
	firstDerivative(mtx, dst, tmp);
}


//second order derivative (central differential quotient) for matrices
template <typename T> void NuLLProcessing::secondDerivative(const Matrix<T>& mtx, Matrix<T>& dst)
{
    uint width = mtx.width();
    uint height = mtx.height();
    T gradMag;

	Matrix<T> dxx(width, height);
	Matrix<T> dyy(width, height);

	Vector<T> stencil(3);
	stencil[0] = 1;
	stencil[1] = -2;
	stencil[2] = 1;

	convolveX(mtx, dxx, stencil);
	convolveY(mtx, dyy, stencil);

    for(uint y=0; y<height; ++y)
    {
        for(uint x=0; x<width; ++x)
        {
            gradMag = sqrt(dxx(x,y)*dxx(x,y)+dyy(x,y)*dyy(x,y));
            dst(x,y) = gradMag;
        }
    }
}

template <typename T> void NuLLProcessing::laplacian(const Matrix<T>& mtx, Matrix<T>& dst)
{
    uint  width = mtx.width();
	uint height = mtx.height();

	Matrix<T> dxx(width, height);
	Matrix<T> dyy(width, height);

	Vector<T> stencil(3);
	stencil[0] = 1;
	stencil[1] = -2;
	stencil[2] = 1;

	convolveX(mtx, dxx, stencil);
	convolveY(mtx, dyy, stencil);

    for(uint y=0; y<height; ++y)
        for(uint x=0; x<width; ++x)
            dst(x,y) = dxx(x,y) + dyy(x,y);
}

//quite similar to an edge detector
//however, bigger neighborhoods can be respected
template <typename T> void NuLLProcessing::localVariance(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	const int radsq = radius * radius;
	int width = mtx.width();
    int height = mtx.height();
	T mean, var;

	std::vector<T> neighbors;
	neighbors.reserve((radius+1)*(radius+1));

	for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if(kx*kx + ky*ky > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }

			mean = 0;
			for(uint i=0; i<neighbors.size(); ++i)
				mean += neighbors[i];

			mean /= neighbors.size();

			var = 0;
			for(uint i=0; i<neighbors.size(); ++i)
				var += (neighbors[i] - mean) * (neighbors[i] - mean);

			var /= neighbors.size();
			dst(x,y) = sqrt(var);

            neighbors.clear();
        }
    }
}

//curvature filter
//useful for corner detection
template <typename T> void NuLLProcessing::meanCurvature(const Matrix<T>& mtx, Matrix<T>& dst)
{
	uint width = mtx.width();
	uint height = mtx.height();
	T dxSq, dySq;
	T k1, k2;

	Matrix<T> dx(width, height);
	Matrix<T> dy(width, height);

	Vector<T> stencil(3);
	stencil[0] = -.5;
	stencil[1] = 0;
	stencil[2] = .5;

	convolveX(mtx, dx, stencil);
	convolveY(mtx, dy, stencil);

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			dxSq = dx(x,y)*dx(x,y);
			dySq = dy(x,y)*dy(x,y);
			k1 = .5*(dxSq + dySq) + sqrt(((dxSq+dySq)/4) - dxSq*dySq + 2*dx(x,y)*dy(x,y));
			k2 = .5*(dxSq + dySq) - sqrt(((dxSq+dySq)/4) - dxSq*dySq + 2*dx(x,y)*dy(x,y));
			dst(x,y) = (k1+k2)/2;
		}
	}
}

template <typename T> void NuLLProcessing::harrisCorners(const Matrix<T>& mtx, Matrix<T>& dst)
{
	uint width = mtx.width();
	uint height = mtx.height();
	T det, tr;

	Matrix<T> dx(width, height);
	Matrix<T> dy(width, height);

	Vector<T> stencil(3);
	stencil[0] = -.5;
	stencil[1] = 0;
	stencil[2] = .5;

	convolveX(mtx, dx, stencil);
	convolveY(mtx, dy, stencil);

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			det = dx(x,y)*dx(x,y) * dy(x,y)*dy(x,y) - 2*dx(x,y)*dy(x,y);
			tr = dx(x,y)*dx(x,y) + dy(x,y)*dy(x,y);
			dst(x,y) = (tr==0)? det : (det/tr);
		}
	}
}

//fixes values which are below/ over given range
template <typename T> void NuLLProcessing::snap(const Matrix<T>& mtx, Matrix<T>& dst, T min, T max)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = std::min<T>(std::max<T>(mtx(x,y), min), max);
}

//filters out entries below threshold
template <typename T> void NuLLProcessing::thresholding(const Matrix<T>& mtx, Matrix<T>& dst, double threshold, double gmin, double gmax)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = (mtx(x,y) >= threshold)? gmax : gmin;
}

//filter entries outside range
template <typename T> void NuLLProcessing::doubleThresholding(const Matrix<T>& mtx, Matrix<T>& dst, double thresholdLower, double thresholdUpper, double gmin, double gmax)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = (mtx(x,y) >= thresholdLower && mtx(x,y) <= thresholdUpper)? gmax : gmin;
}

//automatically calculate threshold value and filter out entries
//uses otsu's method
template <typename T> void NuLLProcessing::automatedThresholding(const Matrix<T>& mtx, Matrix<T>& dst, double gmin, double gmax)
{
	double threshold = 0;

	double curVar = 0;
	double bestVar = 0;

	double cumulMean = 0;
	double totalMean = 0;

	double cumulProb = 0;
	Vector<double> prob(256);

	//probabilities for each colour
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			++prob[(int)mtx(x,y)];

	prob /= (double)(mtx.height() * mtx.width());

	//total mean calculation
	for(uint i=0; i<prob.size(); ++i)
		totalMean += (i+1) * prob[i];

	//find best threshold
	for(double t=0; t<255; ++t)
	{
		cumulProb += prob[(uint)t];
		cumulMean += (t+1) * prob[(int)t];

		curVar = totalMean * cumulProb - cumulMean;
		curVar *= totalMean * cumulProb - cumulMean;
		curVar /= cumulProb * (1 - cumulProb);
		curVar = (cumulProb==0)? 0 : curVar;		//fix division by zero

		if(curVar > bestVar)
		{
			bestVar = curVar;
			threshold = t;
		}
	}

	NuLLProcessing::thresholding(mtx, dst, threshold, gmin, gmax);
}

//simple gamma correction
template <typename T> void NuLLProcessing::gammaCorrection(const Matrix<T>& mtx, Matrix<T>& dst, double gamma)
{
	T maxVal = NuLLTools::maxValue(mtx);
	T gammaInv = 1.0 / gamma;

	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = maxVal * pow((mtx(x,y) / maxVal), gammaInv);
}

//contrast enhancement
template <typename T> void NuLLProcessing::logDynamicCompression(const Matrix<T>& mtx, Matrix<T>& dst, double c)
{
	if(!c)
	{
		c = 255;
		c /= log(1 + NuLLTools::maxValue(mtx));
	}

	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = c * log(1 + mtx(x,y));
}

//rescales all values to fit given range
template <typename T> void NuLLProcessing::affineRescale(const Matrix<T>& mtx, Matrix<T>& dst, double minVal, double maxVal)
{
	T val;
	T imgMax = NuLLTools::maxValue(mtx);
	T imgMin = NuLLTools::minValue(mtx);
	T imgDiff = (imgMax == imgMin)? imgMin : (imgMax - imgMin);		//if imgMax==imgMin, we get div by zero
	T valDiff = maxVal - minVal;

	for(uint y=0; y<mtx.height(); ++y)
	{
		for(uint x=0; x<mtx.width(); ++x)
		{
			val = (mtx(x,y) - imgMin) * (valDiff);
			val /= imgDiff;
			val += minVal;
			dst(x,y) = val;
		}
	}
}

//simple affine transformation
template <typename T> void NuLLProcessing::affineTransform(const Matrix<T>& mtx, Matrix<T>& dst, double a, double b)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = a * mtx(x,y) + b;
}

//downsampling of image
//low-pass filtering before use is recommended
template <typename T> void NuLLProcessing::downsample(const Matrix<T>& mtx, Matrix<T>& dst, uint factorX, uint factorY)
{
	//only one factor given: keep ratio
	if(!factorY)
		factorY = factorX;

	T val;
	const uint width = mtx.width() / factorX;
	const uint height = mtx.height() / factorY;
	const uint scale = factorX * factorY;
	dst.resize(width, height);

	for(uint y=0; y<height; y++)
	{
		for(uint x=0; x<width; x++)
		{
			val = 0;
			for(uint fy=0; fy<factorY; ++fy)
			{
				for(uint fx=0; fx<factorX; ++fx)
				{
					val += mtx(factorX*x + fx, factorY*y + fy);
				}
			}
			dst(x,y) = val / (T)scale;
		}
	}
}

//upscaling of image
template <typename T> void NuLLProcessing::upsample(const Matrix<T>& mtx, Matrix<T>& dst, uint factorX, uint factorY)
{
	if(!factorY)
		factorY = factorX;

	const uint width = mtx.width();
	const uint height = mtx.height();
	dst.resize(width*factorX, height*factorY);

	for(uint y=0; y<height; y++)
	{
		for(uint x=0; x<width; x++)
		{
			for(uint fy=0; fy<factorY; ++fy)
			{
				for(uint fx=0; fx<factorX; ++fx)
				{
					dst(factorX*x + fx, factorY*y + fy) = mtx(x,y);
				}
			}
		}
	}
}

//mirror along x axis
template <typename T> void NuLLProcessing::mirrorX(const Matrix<T>& mtx, Matrix<T>& dst)
{
	const uint width = mtx.width();
	const uint height = mtx.height();

	for(uint y=0; y<height; ++y)
		for(uint x=0; x<width; ++x)
			dst(x,y) = mtx(width-x-1, y);
}

//mirror along y axis
template <typename T> void NuLLProcessing::mirrorY(const Matrix<T>& mtx, Matrix<T>& dst)
{
	const uint width = mtx.width();
	const uint height = mtx.height();

	for(uint y=0; y<height; ++y)
		for(uint x=0; x<width; ++x)
			dst(x,y) = mtx(x, height-y-1);
}

template <typename T> void NuLLProcessing::localMaxima(const Matrix<T>& mtx, Matrix<T>& dst, const int radius)
{
	const uint width = mtx.width();
	const uint height = mtx.height();

	const int radSq = radius * radius;

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			dst(x,y) = mtx(x,y);
			for(int ky=-radius; ky<=radius; ++ky)
			{
				for(int kx=-radius; kx<=radius; ++kx)
				{
					//disc-shaped neighborhood
					if((kx*kx + ky*ky) > radSq)
						continue;

					if(mtx(x,y)<mtx.getMirrored(x+kx,y+ky))
					{
						dst(x,y) = 0;
						break;
					}
				}
			}
		}
	}
}

#endif

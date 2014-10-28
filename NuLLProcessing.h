#ifndef NULLPROCESSING_H
#define NULLPROCESSING_H

#include <vector>
#include <math.h>

#include "Matrix.h"

#include "NuLLTools.inl"
#include "NuLLConvolve.inl"

/*
 *
 * methods for data and signal processing
 * includes:    blurring
 *				point operations
 *              derivative filters (first order, second order)
 *				scaling
 *
 */

namespace NuLLProcessing
{
	//simple blur techniques
	template <typename T> void boxBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void gaussianBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1, double sigma = 0.0);

	template <typename T> void differenceOfGaussians(const Matrix<T>& mtx, Matrix<T>& dst, int sigma1 = 0.0, int sigma2 = 10.0);

	//derivative filters
	template <typename T> void firstDerivative(const Matrix<T>& mtx, Matrix<T>& dst);
	template <typename T> void firstDerivative(const Matrix<T>& mtx, Matrix<T>& dstGrad, Matrix<T>& dstDir);
	template <typename T> void secondDerivative(const Matrix<T>& mtx, Matrix<T>& dst);
	template <typename T> void laplacian(const Matrix<T>& mtx, Matrix<T>& dst);

	//others
	template <typename T> void localVariance(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void meanCurvature(const Matrix<T>& mtx, Matrix<T>& dst);
	template <typename T> void harrisCorners(const Matrix<T>& mtx, Matrix<T>& dst);

	//point operations
	template <typename T> void gammaCorrection(const Matrix<T>& mtx, Matrix<T>& dst, double gamma = 1.2);
	template <typename T> void logDynamicCompression(const Matrix<T>& mtx, Matrix<T>& dst, double c = 0);
	template <typename T> void affineRescale(const Matrix<T>& mtx, Matrix<T>& dst, double minVal = 0, double maxVal = 255);
	template <typename T> void affineTransform(const Matrix<T>& mtx, Matrix<T>& dst, double a = 1, double b = 0);

	template <typename T> void snap(const Matrix<T>& mtx, Matrix<T>& dst, T min = 0, T max = 255);
	template <typename T> void thresholding(const Matrix<T>& mtx, Matrix<T>& dst, double threshold, double gmin = 0, double gmax = 255);
	template <typename T> void doubleThresholding(const Matrix<T>& mtx, Matrix<T>& dst, double thresholdLower, double thresholdUpper, double gmin = 0, double gmax = 255);
	template <typename T> void automatedThresholding(const Matrix<T>& mtx, Matrix<T>& dst, double gmin = 0, double gmax = 255);

	//up- and downscaling
	template <typename T> void downsample(const Matrix<T>& mtx, Matrix<T>& dst, uint factorX=2, uint factorY=0);
	template <typename T> void upsample(const Matrix<T>& mtx, Matrix<T>& dst, uint factor=2, uint factorY=0);

	//mirroring
	template <typename T> void mirrorX(const Matrix<T>& mtx, Matrix<T>& dst);
	template <typename T> void mirrorY(const Matrix<T>& mtx, Matrix<T>& dst);

	//misc
	template <typename T> void localMaxima(const Matrix<T>& mtx, Matrix<T>& dst, const int radius=1);
}

#endif // NULLPROCESSING_H

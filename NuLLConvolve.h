#ifndef NULLCONVOLVE_H
#define NULLCONVOLVE_H

#include "Matrix.h"
#include "Vector.h"

/*
 * Convolution methods for a variety of kernel types
 * Includes kernel generators
 *
 */

namespace NuLLConvolve
{
	//kernel operations
    template <typename T> void normalize(Matrix<T>& mtx);
    template <typename T> void convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Matrix<T>& kernel);
	template <typename T> void convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernel);
	template <typename T> void convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelX, const Vector<T>& kernelY);
	template <typename T> void convolveX(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelX);
	template <typename T> void convolveY(const Matrix<T>& mtx, Matrix<T>& dst, const Vector<T>& kernelY);

	//kernel generators
	template <typename T> void identityKernel(Matrix<T>& dst, const int radius = 1);
	template <typename T> void pillboxKernel(Matrix<T>& dst, const int radius = 1);
    template <typename T> void gaussianKernel(Matrix<T>& dst, const int radius = 1, const double sigma = 0.5);
}

#endif // NULLCONVOLVE_H

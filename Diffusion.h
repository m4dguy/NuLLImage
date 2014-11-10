#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "Matrix.h"
#include "NuLLProcessing.h"

/*
 * Iterative filters using the variational framework
 *
 */

namespace Diffusion
{
	template <typename T> void linearDiffusion(const Matrix<T>& src, Matrix<T>& dst, const double stepsize=.1, const uint iterations=25);
	template <typename T> void nonlinearDiffusion(const Matrix<T>& src, Matrix<T>& dst, const double lambda = 2.5, const double stepsize = .1, const uint iterations = 25);
	template <typename T> void regularizedNonlinearDiffusion(const Matrix<T>& src, Matrix<T>& dst, const double lambda = 2.5, const double sigma = 1., const double stepsize = .1, const uint iterations = 25);
}

#endif
#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "Matrix.h"

/*
 * Iterative filters using the variational framework
 *
 */

namespace Diffusion
{
	template <typename T> void linearDiffusion(const Matrix<T>& src, Matrix<T>& dst, const double stepsize=.1, const uint iterations=25);
}

#endif
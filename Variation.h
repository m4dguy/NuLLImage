#ifndef VARIATION_H
#define VARIATION_H

#include "Matrix.h"

/*
 * Iterative filters using the variational framework
 *
 */

namespace Variation
{
	template <typename T> void variationalBlur(const Matrix<T>& src, Matrix<T>& dst, const double alpha=1., const uint iterations=100);
	//template <typename T> void charbonnierDenoise(const Matrix<T>& src, Matrix<T>& dst, const double alpha=1., const double lambda=1., const uint iterations=100);
}

#endif
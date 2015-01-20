#ifndef FLOW_H
#define FLOW_H

#include "Matrix.h"

/*
 * Optical flow methods
 *
 */

namespace Flow
{
	template <typename T> void hornSchunck(const Matrix<T>& src1, const Matrix<T>& src2, Matrix<T>& dstX, Matrix<T>& dstY, const double alpha = 10., const uint iterations=25);
}

#endif
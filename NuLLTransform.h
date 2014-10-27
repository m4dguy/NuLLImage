#ifndef NULLTRANSFORM_H
#define NULLTRANSFORM_H

#include <math.h>

#include "Utils.h"
#include "Vector.h"
#include "Matrix.h"

namespace NuLLTransform
{
    const double pi = 2.0 * asinf(1.0);

    //standard slow fourier
	template <typename T> void fourierTransform(const Matrix<T>& src, Matrix<T>& dstReal, Matrix<T> dstImag);
	template <typename T> void cosineTransform(const Matrix<T>& src, Matrix<T>& dst);
}

#endif // NULLTRANSFORM_H

#ifndef FLOW_INL
#define FLOW_INL

#include "Flow.h"

//simplified horn schunck optic flow calculation
template <typename T> void Flow::hornSchunck(const Matrix<T>& src1, const Matrix<T>& src2, Matrix<T>& dstX, Matrix<T>& dstY, const double alpha, const uint iterations)
{
	T laplaceu, laplacev;
	const size_t width = src1.width();
	const size_t height = src1.height();

	Matrix<T> dx(width, height);
	Matrix<T> dy(width, height);
	Matrix<T> dt(width, height);

	Matrix<T> tmpX(width, height);
	Matrix<T> tmpY(width, height);

	//gradient calculation
	for (uint y = 0; y < height; ++y)
	{
		for (uint x = 0; x < width; ++x)
		{
			dx(x, y) = .5 * (src1.getMirrored(x + 1, y) - src1.getMirrored(x - 1, y) + src2.getMirrored(x + 1, y) - src2.getMirrored(x - 1, y));
			dy(x, y) = .5 * (src1.getMirrored(x, y + 1) - src1.getMirrored(x, y - 1) + src2.getMirrored(x, y + 1) - src2.getMirrored(x, y - 1));
			dt(x, y) = src2(x,y) - src1(x,y);
		}
	}

	//flow calculation
	for (uint i = 0; i<iterations; ++i)
	{
		for (uint y = 0; y<height; ++y)
		{
			for (uint x = 0; x<width; ++x)
			{
				laplaceu = (dstX.getMirrored(x + 1, y) - 2 * dstX(x, y) + dstX.getMirrored(x - 1, y)) + (dstX.getMirrored(x, y + 1) - 2 * dstX(x, y) + dstX.getMirrored(x, y - 1));
				laplacev = (dstY.getMirrored(x + 1, y) - 2 * dstY(x, y) + dstY.getMirrored(x - 1, y)) + (dstY.getMirrored(x, y + 1) - 2 * dstY(x, y) + dstY.getMirrored(x, y - 1));

				tmpX(x, y) = (-dx(x, y)*dt(x, y) - (dx(x, y)*dy(x, y)*dstY(x, y) - alpha*laplaceu)) / (dx(x, y)*dx(x, y) + alpha*8.);
				tmpY(x, y) = (-dy(x, y)*dt(x, y) - (dx(x, y)*dy(x, y)*dstX(x, y) - alpha*laplacev)) / (dy(x, y)*dy(x, y) + alpha*8.);
			}
		}
		dstX.swap(tmpX);
		dstY.swap(tmpY);
	}

}


#endif
#ifndef NULLMORPHOLOGY_H
#define NULLMORPHOLOGY_H

#include <algorithm>
#include <vector>
#include <math.h>

/*
 *
 *
 */

namespace NuLLMorphology
{
	//simple morphological filters
    template <typename T> void medianFilter(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1, float percentile = .5f);
    template <typename T> void dilation(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
    template <typename T> void erosion(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
    template <typename T> void opening(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
    template <typename T> void closing(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void whiteTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void blackTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void selfdualTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);

	//non-flat morphology
	template <typename T> void distanceTransform(const Matrix<T>& mtx, Matrix<T>& dst);
}

#include NuLLMorphology.inl"

#endif // NULLPROCESSING_H

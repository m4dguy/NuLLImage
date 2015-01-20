#ifndef NULLFILE_H
#define NULLFILE_H

#include <fstream>

#include "Matrix.h"

namespace NuLLFile
{
	template<typename T> int writeCSV(const Matrix<T>& src, const char* filename);

	template<typename T> int writePGM(const Matrix<T>& src, const char* filename);
	template<typename T> int writePPM(const Matrix<T>& srcR, const Matrix<T>& srcG, const Matrix<T>& srcB, const char* filename);
};

#include "NuLLFile.inl"

#endif // NULLFILE_H

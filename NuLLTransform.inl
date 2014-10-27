#ifndef NULLTRANSFORM_INL
#define NULLTRANSFORM_INL

#include "NuLLTransform.h"

//can be used for matrix conversion
//copy content of matrix src to matrix dst
//useful for type conversion
template <typename T> void NuLLTransform::fourierTransform(const Matrix<T>& src, Matrix<T>& dstReal, Matrix<T> dstImag)
{
    const size_t width = src.width();
    const size_t height = src.height();

    T realPart, imaginaryPart;
    T sumMr, sumNr;
    T sumMi, sumNi;

    const T W = src.width();                    //avoid casting
    const T H = src.height();                   //avoid casting
    const T scaling = 1./sqrt(W * H);
    const T pi2 = 2. * pi;

    for(uint y=0; y<height; ++y)
    {
        for(uint x=0; x<width; ++x)
        {
            sumMr = 0;
            sumMi = 0;
            for(uint m=0; m<height; ++m)
            {
                sumNr = 0;
                sumNi = 0;
                for(uint n=0; n<width; ++n)
                {
                    //sumN += src(x,y) * exp((pi2*x*m)/h) * exp(-(pi2*y*n)/w);
                    sumNr += src(x,y) * (cos((pi2*x*m)/H) * cos((pi2*y*n)/W) - sin((pi2*x*m)/H) * sin((pi2*y*n)/W));
                    sumNi += src(x,y) * (sin((pi2*x*m)/H) * cos((pi2*y*n)/W) + sin((pi2*y*n)/W) * cos((pi2*x*m)/H)) * (-1);
                }
                sumMr += sumNr;
                sumMi += sumNi;
            }
            dstReal(x,y) = scaling * sumMr;
            dstImag(x,y) = scaling * sumMi;
        }
    }
}

template <typename T> void NuLLTransform::cosineTransform(const Matrix<T>& src, Matrix<T>& dst)
{
    const size_t width = src.width();
    const size_t height = src.height();

    T sum;

    const T scaleX0 = sqrt(1. / width);
    const T scaleY0 = sqrt(1. / height);

    const T scaleXn = sqrt(2. / width);
    const T scaleYn = sqrt(2. / height);

    const T w2 = 2 * width;
    const T h2 = 2 * height;

    T scaleX, scaleY;

    for(uint y=0; y<height; ++y)
    {
        scaleY = y? scaleYn : scaleY0;
        for(uint x=0; x<width; ++x)
        {
            scaleX = x? scaleXn : scaleX0;
            sum = 0;
            for(uint m=1; m<height; ++m)
            {
                for(uint n=1; n<width; ++n)
                {
                    sum += src(m,n) * scaleX * scaleY * cos((pi*(2*m+1)*x) / w2) * cos((pi*(2*n+1)*y) / h2);
                }
            }
            dst(x,y) = sum;
        }
    }
}

#endif

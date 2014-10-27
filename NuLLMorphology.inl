#ifndef NULLMORPHOLOGY_INL
#define NULLMORPHOLOGY_INL

#include "NuLLMorphology.h"

//median filtering for matrices
template <typename T> void NuLLMorphology::medianFilter(const Matrix<T>& mtx, Matrix<T>& dst, int radius, float percentile)
{
    const int radsq = radius * radius;
    int width = mtx.width();
    int height = mtx.height();

    std::vector<T> neighbors;
	neighbors.reserve((radius+1)*(radius+1));

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if((kx*kx + ky*ky) > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }
            std::sort(neighbors.begin(), neighbors.end());
            dst(x,y) = neighbors[(uint)((neighbors.size()-1) * percentile)];
            neighbors.clear();
        }
    }
}

//dilation for matrices
template <typename T> void NuLLMorphology::dilation(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    int width = mtx.width();
    int height = mtx.height();

    T sup;
    std::vector<T> neighbors;
    neighbors.reserve((radius+1)*(radius+1));
    const int radsq = radius * radius;

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if(kx*kx + ky*ky > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }

            sup = neighbors[0];
            for(uint i=1; i<neighbors.size(); ++i)
            {
                if(neighbors[i] < sup)
                    sup = neighbors[i];
            }
            dst(x,y) = sup;
            neighbors.clear();
        }
    }
}

//erosion for matrices
template <typename T> void NuLLMorphology::erosion(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    int width = mtx.width();
    int height = mtx.height();

    T inf;
    std::vector<T> neighbors;
    neighbors.reserve((radius+1)*(radius+1));
    const int radsq = radius * radius;

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if(kx*kx + ky*ky > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }

            inf = neighbors[0];
            for(uint i=1; i<neighbors.size(); ++i)
            {
                if(neighbors[i] > inf)
                    inf = neighbors[i];
            }
            dst(x,y) = inf;
            neighbors.clear();
        }
    }
}

//opening for matrices
template <typename T> void NuLLMorphology::opening(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    Matrix<T> tmp(mtx.width(), mtx.height());
    erosion(mtx, tmp, radius);
    dilation(tmp, dst, radius);
}

//closing for matrices
template <typename T> void NuLLMorphology::closing(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    Matrix<T> tmp(mtx.width(), mtx.height());
    dilation(mtx, tmp, radius);
    erosion(tmp, dst, radius);
}

//white tophat operation for matrices
template <typename T> void NuLLMorphology::whiteTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
    opening(mtx, tmp, radius);
    dst = mtx;
	dst -= tmp;
}

//black tophat operation for matrices
template <typename T> void NuLLMorphology::blackTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
	closing(mtx, dst, radius);
	dst -= mtx;
}

template <typename T> void NuLLMorphology::selfdualTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
	blackTopHat(mtx, dst, radius);
	whiteTopHat(mtx, tmp, radius);
	dst += tmp;
}

//distance transform
//uses fast scheme vd boomgaard
//use on binary images with values: 0, 255
template <typename T> void NuLLMorphology::distanceTransform(const Matrix<T>& mtx, Matrix<T>& dst)
{
	uint  width = mtx.width();
	uint height = mtx.height();

	//max distance
	const T inf = (width * width + height * height);

	uint in;						//point with distance information
	T val, valMin;					//calculated distance

	Matrix<T> tmp(width, height);
	Matrix<T> res(width, height);		//starts with zero distance

	//init distances
	for (uint y = 0; y<height; ++y)
	{
		for (uint x = 0; x<width; ++x)
		{
			if (mtx(x, y))
				tmp(x, y) = 0;
			else
				tmp(x, y) = inf;
		}
	}

	//transform in x direction
	for (uint y = 0; y<height; ++y)
	{
		res(0, y) = inf;

		//forward direction
		for (uint x = 1; x<width; ++x)
		{
			//find reference point inside object
			if (tmp(x, y) <= res(x - 1, y))
			{
				res(x, y) = tmp(x, y);
				in = x;
			}
			else
			{
				valMin = inf;
				for (uint xi = in; xi<x; ++xi)
				{
					val = tmp(xi, y) + (x - xi) * (x - xi);
					if (val <= valMin)
					{
						valMin = val;
						in = xi;
					}
				}
				res(x, y) = valMin;
			}
		}

		//backward direction
		for (int x = width - 2; x >= 0; --x)
		{
			//find reference point inside object
			if (res(x, y) <= tmp(x + 1, y))
			{
				tmp(x, y) = res(x, y);
				in = x;
			}
			else
			{
				valMin = inf;
				for (int xi = in; xi>x; --xi)
				{
					val = res(xi, y) + (x - xi) * (x - xi);
					if (val <= valMin)
					{
						valMin = val;
						in = xi;
					}
				}
				tmp(x, y) = valMin;
			}
		}
	}

	//transform in y direction
	for (uint x = 0; x<width; ++x)
	{
		//forward direction
		for (uint y = 1; y<height; ++y)
		{
			if (tmp(x, y) <= res(x, y - 1))
			{
				res(x, y) = tmp(x, y);
				in = y;
			}
			else
			{
				valMin = inf;
				for (uint yi = in; yi<y; ++yi)
				{
					val = tmp(x, yi) + (y - yi) * (y - yi);
					if (val <= valMin)
					{
						valMin = val;
						in = yi;
					}
				}
				res(x, y) = valMin;
			}
		}

		//backward direction
		for (int y = height - 2; y >= 0; --y)
		{
			if (res(x, y) <= tmp(x, y + 1))
			{
				tmp(x, y) = res(x, y);
				in = y;
			}
			else
			{
				valMin = inf;
				for (int yi = in; yi>y; --yi)
				{
					val = res(x, yi) + (y - yi) * (y - yi);
					if (val <= valMin)
					{
						valMin = val;
						in = yi;
					}
				}
				tmp(x, y) = valMin;
			}
		}
	}
	dst.swap(tmp);
}

#endif

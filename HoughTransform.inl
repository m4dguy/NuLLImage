#include "HoughTransform.h"

//hough transform for circles
//make sure img contains edge information
template<typename T> void HoughTransform::houghCircles(const Matrix<T>& img, std::vector<HoughCircle>& circles, double threshold, int minRad, int maxRad)
{
	HoughCircle circ;
	int rx, ry;
	uint pixels;
	double score;
	double radSq, radSqLow, radSqUp;

	const uint width = img.width();
	const uint height = img.height();

	std::vector<std::pair<uint,uint> > mask;
	mask.reserve(maxRad * maxRad);			//TODO: too many reserved; #pixels <= ceil(#area / 2)

	for(int rad=minRad; rad<=maxRad; ++rad)
	{
		radSqLow = (rad - .5) * (rad - .5);
		radSqUp = (rad + .5) * (rad + .5);

		//precalculate pixels in circle
		for(rx=-rad; rx<=rad; ++rx)
		{
			for(ry=-rad; ry<=rad; ++ry)
			{
				radSq = (rx*rx + ry*ry);
				if((radSqLow <= radSq) && (radSq < radSqUp))
				{
					mask.push_back(std::make_pair(rx,ry));
				}
			}
		}
		pixels = mask.size();

		//check image
		for(uint y=rad; y<height-rad; ++y)
		{
			for(uint x=rad; x<width-rad; ++x)
			{
				for(uint i=0; i<pixels; ++i)
				{
					rx = mask[i].first;
					ry = mask[i].second;
					score += img(x+rx, y+ry);
				}
				
				//check if circle found
				score /= pixels;
				if(score >= threshold)
				{
					circ.x = x;
					circ.y = y;
					circ.radius = rad;
					circ.score = score;
					circles.push_back(circ);
				}
			}
		}
		mask.clear();
	}
}

template<typename T> void HoughTransform::plotCircles(const std::vector<HoughCircle>& circles, Matrix<T>& img, uint radius)
{
	uint x, y;
	int rad, radSq;
	double radSqLow, radSqUp;

	const uint width = img.width();
	const uint height = img.height();

	for(uint i=0; i<circles.size(); ++i)
	{
		for(uint r=0; r<radius; ++r)
		{
			HoughCircle c = circles[i];
			x = c.x;
			y = c.y;
			rad = c.radius - r;
			radSqLow = (rad - .5) * (rad - .5);
			radSqUp = (rad + .5) * (rad + .5);

			for(int ry=-rad; ry<=rad; ++ry)
			{
				for(int rx=-rad; rx<=rad; ++rx)
				{
					radSq = (rx*rx + ry*ry);
					if((radSqLow <= radSq) && (radSq < radSqUp))
					{
						img(x+rx, y+ry) = 255;
					}
					
				}
			}
		}
	}
}

template <typename T> inline void HoughTransform::plotCirclesFilled(const std::vector<HoughCircle>& circles, Matrix<T>& img)
{
	uint x, y;
	int rad, radSq;
	const uint width = img.width();
	const uint height = img.height();

	for(uint i=0; i<circles.size(); ++i)
	{
		HoughCircle c = circles[i];
		x = c.x;
		y = c.y;
		rad = c.radius;
		radSq = rad * rad;

		for(int ry=-rad; ry<=rad; ++ry)
		{
			for(int rx=-rad; rx<=rad; ++rx)
			{
				if((rx*rx + ry*ry) <= radSq)
				{
					img(x+rx, y+ry) = circles[i].marker;
				}
			}
		}
	}
}

//merge nearby circles
//the circle with the higher hough score is kept
void HoughTransform::mergeCircles(std::vector<HoughCircle>& circles, uint maxDist)
{
	uint dist;
	uint bSize;
	const uint maxDistSq = maxDist * maxDist;

	HoughCircle c;

	std::vector<uint> done;
	std::vector<HoughCircle> batch;
	std::vector<HoughCircle> res;
	done.resize(circles.size());
	batch.reserve(circles.size());
	res.reserve(circles.size());

	for(uint i=0; i<circles.size()-1; ++i)
	{
		if(done[i])
			continue;

		for(uint j=i+1; j<circles.size(); ++j)
		{
			if(done[j])
				continue;

			dist = ((circles[i].x - circles[j].x) * (circles[i].x - circles[j].x)) + ((circles[i].y - circles[j].y) * (circles[i].y - circles[j].y));
			if(dist <= maxDistSq)
			{
				done[i] = 1;
				done[j] = 1;
				batch.push_back(circles[j]);
			}
		}

		bSize = batch.size();
		if(bSize)
		{
			c = batch[0];
			for(uint b=0; b<bSize; ++b)
			{
				if(batch[b].score > c.score)
					c = batch[b];
			}
			res.push_back(c);
			batch.clear();
		}
	}
	res.swap(circles);
}

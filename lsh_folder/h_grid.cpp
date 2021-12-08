//file:h_grid.cpp//
#include "h_grid.hpp"
#include "params.hpp"
#include "object.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <vector>

h_grid::h_grid()
{	
	// picks single precision real t uniformly in [0,delta)^2
	const float lower_bound = 0.0;
	const float upper_bound = (float) delta;
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> distr(lower_bound, upper_bound);
	
	t1 = distr(generator);
	t2 = distr(generator);
}

Object * h_grid::operator()(const Object& time_series_curve) const
{
	// first of all the function snaps the given Object-time series to grid, without multiplying by delta or shifting by t
	// we do this so that we can remove duplicate points from snapping by comparing integers (otherwise we would compare floats , oof)
	
	// vector holds snapped points-coordinates
	std::vector <int> snapped_curve = time_series_curve.snap();

	// vector of grid curve coordinates after removing duplicates from snapping
	std::vector <float> grid_curve = time_series_curve.remove_dupls(snapped_curve);

	// multiply by delta and shift by t to get final grid curve
	if (metric_func == "discrete")		// 2d grid
	{
		for (int i = 0; i < (int) grid_curve.size(); i+=2)
		{
			grid_curve[i] = grid_curve[i] * delta + t1;
			grid_curve[i+1] = grid_curve[i+1] * delta + t2;
		}
	}
	else if (metric_func == "continuous")	// 1d grid
	{
		for (int i = 0; i < (int) grid_curve.size(); ++i)
			grid_curve[i] = grid_curve[i] * delta + t1;
	}
	
	// do padding necessary
	time_series_curve.pad(grid_curve);

	// return the new time series
	if (metric_func == "discrete")
		return new time_series(grid_curve);
	else // (metric == "continuous")
		return new Object(grid_curve);
}
	

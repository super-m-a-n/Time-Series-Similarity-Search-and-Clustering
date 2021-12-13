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
	t.push_back(distr(generator));
	t.push_back(distr(generator));

}

Abstract_Object * h_grid::operator()(const Abstract_Object& time_series) const
{
	return time_series.to_grid_curve(this->t);
}
	

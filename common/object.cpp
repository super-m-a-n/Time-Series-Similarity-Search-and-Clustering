//file:object.cpp//
#include "object.hpp"
#include "params.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <vector>



/////////////////////////////// CLASS OBJECT ///////////////////////////////////////////

///////////////////////// CONSTRUCTION / DESTRUCTION ///////////////////////////////////

Object::Object()
{
	float norm_squared = 0.0;
	int dim;
	
	// different Object size depending on algorithm
	if (algorithm == "LSH" || algorithm == "Hypercube")
		dim = d;
	else if (algorithm == "Frechet" && metric_func == "discrete")
		dim = 2*d;
	else if (algorithm == "Frechet" && metric_func == "continuous")
		dim = d;

	for (int i = 0; i < dim; ++i)
	{
		std::random_device rd;
		std::mt19937 generator(rd());
		std::normal_distribution<float> distr(0, 1);		// each coordinate follows normal(0,1) distribution

		vector.push_back(distr(generator));

		norm_squared += vector[i] * vector[i];
	}

	for (int i = 0; i < dim; ++i)
		vector[i] = vector[i] / sqrt(norm_squared);			// divide each coordinate by norm, to normalize point-object
}

Object::Object(std::vector <float> & input_vector, std::string & object_name) : identifier(object_name)
{
	for (int i = 0; i < (int) input_vector.size(); ++i)
		vector.push_back(input_vector[i]);
}

Object::Object(std::vector <float> & input_vector)
{
	for (int i = 0; i < (int) input_vector.size(); ++i)
		vector.push_back(input_vector[i]);
}

Object::Object(std::string & object_name) : identifier(object_name){}

Object::~Object() {}

///////////////////////// GETTERS /////////////////////////////////////////////

const std::string & Object::get_name() const
{
	return this->identifier;
}

float Object::get_ith(int i) const
{
	return this->vector[i];
}

int Object::get_dim() const
{
	return (int) (this->vector).size();
}

///////////////////////// SETTERS /////////////////////////////////////////////

void Object::set(const Object& p)
{
	// copy identifier
	this->identifier = p.identifier;

	if (this->get_dim() != p.get_dim())		// object dimensions should match
	{
		std::cerr << "Warning : Object::set : dimensions don't match\n\n";
		return;
	}

	// copy coordinates
	for (int i = 0; i < this->get_dim(); ++i)
		this->vector[i] = p.vector[i];
}


void Object::set_ith(int i, float value)
{
	this->vector[i] = value;
}

////////////////////////  PRINTS ///////////////////////////////////////////////

void Object::print() const
{
	std::cout << "Object " << this->identifier << "--> ";
	
	this->print_coordinates();
	
}

void Object::print_coordinates() const{
	std::cout << "(";
	for (int i = 0; i < (this->get_dim() - 1); ++i)
		std::cout << this->vector[i] << ",";
	std::cout << this->vector[this->get_dim() - 1];
	std::cout << ")";
}

void Object::print(std::ofstream & file) const
{
	file << "Object " << this->identifier << "--> ";
	
	this->print_coordinates(file);
}

void Object::print_coordinates(std::ofstream & file) const{
	file << "(";
	for (int i = 0; i < (this->get_dim() - 1); ++i)
		file << this->vector[i] << ",";
	file << this->vector[this->get_dim() - 1];
	file << ")";
}

/////////////////////// OBJECT OPERATIONS ///////////////////////////////////////////

float Object::inner_prod(const Object& p) const
{
	float inner_prod = 0.0;

	if (this->get_dim() != p.get_dim())		// object dimensions should match for inner product
	{
		std::cerr << "Warning : Object::inner_prod : dimensions don't match\n\n";
		return inner_prod;		//  returns garbage value
	}

	for (int i = 0; i < this->get_dim(); ++i)
		inner_prod += this->vector[i] * p.vector[i];

	return inner_prod;
}

double Object::euclidean_distance(const Object& p) const
{
	double dist_squared = 0.0;

	if (this->get_dim() != p.get_dim())		// object dimensions should match for euclidean distance
	{
		std::cerr << "Warning : Object::euclidean_distance : dimensions don't match\n\n";
		return dist_squared;		//  returns garbage value
	}

	for (int i = 0; i < this->get_dim(); ++i)
		dist_squared += (double) (this->vector[i] - p.vector[i]) * (double) (this->vector[i] - p.vector[i]);

	return sqrt(dist_squared);
}


std::vector <int> Object::snap() const
{
	std::vector <int> snapped_object;

	for (int i = 0; i < this->get_dim(); ++i)
		snapped_object.push_back(floor(this->vector[i] / delta + 1/2));		// snap each coordinate to new integer coordinate of grid

	return snapped_object;
}

std::vector <float> Object::remove_dupls(std::vector <int> & snapped_curve) const
{
	std::vector <float> grid_curve;
	
	// for each point of curve, one boolean variable stating if it is duplicate or not
	std::vector <bool> duplicate (this->get_complexity(), false);

	// now we mark duplicate points from snapping
	int pos = 0;	// index of last non duplicate point
	for (int i = 1; i < this->get_complexity(); ++i)	// for each point of curve
	{
		// check if point i is the same as last non duplicate point
		
		// if metric is discete time series is 2 dimensional
		if (metric_func == "discrete")
		{
			// check if i-th point is equal to last unique point
			if (snapped_curve[2*i] == snapped_curve[2*pos] && snapped_curve[2*i+1] == snapped_curve[2*pos+1])
				duplicate[i] = true;	// if it same as last non duplicate point, then its a duplicate
			else
				pos = i;		// otherwise it becomes the new last non duplicate point
		}
		// if metric is continuous,  time series is 1 dimensional
		else if (metric_func == "continuous")
		{
			// check if i-th point is equal to last non duplicate point
			if (snapped_curve[i] == snapped_curve[pos])
				duplicate[i] = true;
			else
				pos = i;
		}
	}

	// push back the non duplicate points to grid curve vector
	for (int i = 0; i < this->get_complexity(); ++i)
	{
		if (duplicate[i] == false)
		{
			if (metric_func == "continuous")
				grid_curve.push_back(snapped_curve[i]);
			else if (metric_func == "discrete")
			{
				grid_curve.push_back(snapped_curve[2*i]);
				grid_curve.push_back(snapped_curve[2*i+1]);
			}
		}
	}

	return grid_curve;
}

void Object::pad(std::vector <float> & grid_curve) const
{
	// finally check if complexity/dimension of snapped grid time_series diminished, and apply padding if necessary
	// with special big padding number
	unsigned long int M = 1000000;  /// further testing required
	for (int i = (int) grid_curve.size(); i < this->get_dim(); ++i)
		grid_curve.push_back(M);
}

double euclidean(const Object & p, const Object & q)
{
	return p.euclidean_distance(q);
}


/////////////////////////////// CLASS TIME_SERIES ///////////////////////////////////////////

time_series::time_series(std::vector <float> & input_vector, std::string & curve_name) : Object(curve_name), complexity(input_vector.size())
{
	int dim = (int) input_vector.size();
	float x_val = 1;

	// by default time_series object has 2*dim size (dim for x-values, dim for y-values)
	for (int i = 0; i < 2*dim; ++i)
	{
		if (i % 2 == 0)		// even positions hold x-values
		{
			vector.push_back(x_val);
			x_val++;
		}
		else
			vector.push_back(input_vector[(i-1)/2]); 	// odd positions hold y-values
	}
}

time_series::time_series(std::vector <float> & input_vector) : Object(input_vector) , complexity(input_vector.size()/2){}

int time_series::get_complexity() const
{
	return this->complexity;
}

double time_series::discrete_frechet_distance(const Object & P) const
{
	// a vector of vectors that will serve as the 2D array for dynamic programming
	std::vector <std::vector <double> > OPT(this->complexity, std::vector <double> (P.get_complexity()));
	// initialize first square at (0,0)
	OPT[0][0] = norm(this->vector[0], this->vector[1], P.get_ith(0), P.get_ith(1));
	// initialize first column of array
	for (int i = 1; i < this->complexity; i++)
		OPT[i][0] = std::max(OPT[i-1][0], norm(this->vector[2*i], this->vector[2*i+1], P.get_ith(0), P.get_ith(1)));
	// initialize first row of array
	for (int j = 1; j < P.get_complexity(); j++)
		OPT[0][j] = std::max(OPT[0][j-1], norm(this->vector[0], this->vector[1], P.get_ith(2*j), P.get_ith(2*j+1)));

	// initialize rest of array (i > 0 and j > 0)
	for (int i = 1; i < this->complexity; i++)
		for (int j = 1; j < P.get_complexity(); j++)
			OPT[i][j] = std::max(std::min(OPT[i-1][j], std::min(OPT[i-1][j-1], OPT[i][j-1])), norm(this->vector[2*i], this->vector[2*i+1], P.get_ith(2*j), P.get_ith(2*j+1)));


	return OPT[this->complexity - 1][P.get_complexity() - 1];	// value for frechet distance is at top right corner of array
}


double discrete_frechet(const Object & P, const Object & Q)
{
	return P.discrete_frechet_distance(Q);
}

double norm(float x1, float y1, float x2, float y2)
{
	return sqrt( ((double) x1 - (double) x2)  *  ((double) x1 - (double) x2) + ((double) y1 - (double) y2) * ((double) y1 - (double) y2) );
}
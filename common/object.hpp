//file:object.hpp//
#ifndef _OBJECT_HPP_
#define _OBJECT_HPP_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "params.hpp"

// class Object holds the data of the input objects-points of the dataset
// i.e. the coordinates in a d-vector and the object name, from input file

class Object 
{
protected:
	std::string identifier; 		// the object identifier name, as read from input file
	std::vector <float>  vector;	// the coordinates of the d-dimensional point-object in an array

public:
	///////////////////////// CONSTRUCTION / DESTRUCTION //////////////////////////////////////////////////////////////////
	
	// default constructor creates a random normalized d-dimensional point-object, each coordinate follows normal(0,1) distribution 
	Object();
	// constructor through another d-dimensional input array
	Object(std::vector <float> & input_vector, std::string & object_name);
	// constructor through another d-dimensional input array with empty object name
	Object(std::vector <float> & input_vector);
	// Object constructor that just initializes the object identifier (used by potential subclasses)
	Object(std::string & object_name);
	// destructor
	virtual ~Object();

	///////////////////////// GETTERS /////////////////////////////////////////////
	
	// name identifier getter
	const std::string & get_name() const; 
	// object dimensionality getter
	int get_dim() const; 
	// gets ith coordinate of Object
	float get_ith(int i) const;

	///////////////////////// SETTERS /////////////////////////////////////////////
	
	// sets caller object's info to given arg object's info (does a copy basically)
	void set(const Object& p);
	// sets ith coordinate of Object to value
	void set_ith(int i, float value);

	////////////////////////  PRINTS ///////////////////////////////////////////////
	
	// print method for debugging
	void print() const;
	// write object to file
	void print(std::ofstream & file) const;
	//Same as with the above prints but now print only the coordinates of the object
	void print_coordinates() const;
	// same as above, but print to file
	void print_coordinates(std::ofstream & file) const;

	/////////////////////// OBJECT OPERATIONS ///////////////////////////////////////////
	
	// calculates the inner-product of calling object with given object p (both d-dimensional)
	float inner_prod(const Object& p) const;
	// calculates the euclidean distance between caller object and argument object
	double euclidean_distance(const Object& p) const;
	// computes the discrete frechet distance between caller object and argument object
	virtual double discrete_frechet_distance(const Object & P) const {return 0.0;}		// for plain object frechet_distance is not defined (zero) TODO :maybe change this
	// returns objects complexity
	virtual int get_complexity() const {return this->get_dim();}	// for plain object complexity == dimension
	// returns a vector with snapped coordinates to grid integers
	std::vector <int> snap() const;
	// returns a vector with snapped coordinates to grid integers after removing
	std::vector <float> remove_dupls(std::vector <int> & snapped_curve) const;
	// pads given vector of grid coordinates if necessary
	void pad(std::vector <float> & grid_curve) const;

};


// class for a 2D time_series object.  We treat a 2D time_series object as an Object essentially, with double size by default (both x,y values)
class time_series : public Object
{
private:
	int complexity;		// the complexity (number of points) of time_series
public:
	///////////////////////// CONSTRUCTION / DESTRUCTION ////////////////////////////////

	// constructor through another d-dimensional input array with just y-values
	time_series(std::vector <float> & input_vector, std::string & curve_name);
	// constructor through another input array with both x-values and y-values
	time_series(std::vector <float> & input_vector);
	
	/////////////////////// TIME_SERIES OPERATIONS ///////////////////////////////////////////

	// computes the discrete frechet distance between caller time series and argument time series
	double discrete_frechet_distance(const Object & P) const;
	// returns curve's complexity
	int get_complexity() const;
};

// metric wrappers
double euclidean(const Object & p, const Object & q);
double discrete_frechet(const Object & P, const Object & Q);
double norm(float x1, float y1, float x2, float y2);

#endif
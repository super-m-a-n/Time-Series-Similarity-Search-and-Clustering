#include "acutest.h"
#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <time.h>
#include <cstdlib>
#include "object.hpp"

// global program parameters (needed to be defined because of params.hpp)
int N = 1, R, d = 0, w, n = 0;		// global
int k, L;							// vector lsh
int d1, probes, M;					// vector hypercube
double delta, epsilon = 0.01;		// curve lsh (frechet)
std::string algorithm, metric_func;	// variable for algorithm , metric used for frechet


bool compare_traversals(std::list<std::pair<int, int>> traversal1, std::list<std::pair<int, int>> traversal2){
    if (traversal1.size() != traversal2.size()) return false;
    
    for (int i = 0 ; i < (int) traversal1.size() ; i++){
        std::pair<int,int> a = traversal1.front();
        std::pair<int,int> b = traversal2.front();

        if (a.first != b.first || a.second != b.second) return false;
        traversal1.pop_front();
        traversal2.pop_front();
    }

    return true;
}

bool compare_mean_curves(std::vector<std::pair<float, float>> mean_curve1, std::vector<std::pair<float, float>> mean_curve2){
    if (mean_curve1.size() != mean_curve2.size()) return false;
    
    for (int i = 0 ; i < (int) mean_curve1.size() ; i++){
        std::pair<float,float> a = mean_curve1[i];
        std::pair<float,float> b = mean_curve2[i];

        if (std::abs(a.first - b.first) > 0.00001 || std::abs(a.second - b.second) > 0.00001) return false;
    }

    return true;
}


void testing_search(void){
    std::vector<float> curve1 {1, 4,    10, 1000, 6, 10};
    std::vector<float> curve2 {1, 4, 6, 12, 1000, 6, 11};
    Object object1(curve1);
    Object object2(curve2);

    time_series t_series1 (curve1);
    time_series t_series2 (curve2);

    TEST_CHECK(std::abs(discrete_frechet(object1, object2) - 2) < 0.00001); 
    TEST_CHECK(std::abs(discrete_frechet(t_series1, t_series2) - sqrt(5)) < 0.00001);

    std::list<std::pair<int, int>> best_traversal;

    best_traversal.push_back (std::make_pair(0,0));
    best_traversal.push_back (std::make_pair(1,1));
    best_traversal.push_back (std::make_pair(1,2));
    best_traversal.push_back (std::make_pair(2,3));
    best_traversal.push_back (std::make_pair(3,4));
    best_traversal.push_back (std::make_pair(4,5));
    best_traversal.push_back (std::make_pair(5,6));

    TEST_CHECK(compare_traversals(best_traversal, t_series1.best_traversal(&t_series2)));
    
    std::vector <std::pair <float, float> > mean_curve;
    mean_curve.push_back (std::make_pair(1.0,1.0));
    mean_curve.push_back (std::make_pair(2.0,4.0));
    mean_curve.push_back (std::make_pair(2.5,5.0));
    mean_curve.push_back (std::make_pair(3.5,11.0));
    mean_curve.push_back (std::make_pair(4.5,1000.0));
    mean_curve.push_back (std::make_pair(5.5,6.0));
    mean_curve.push_back (std::make_pair(6.5,10.5));

    TEST_CHECK(compare_mean_curves(mean_curve, t_series1.mean_curve_without_filtering(&t_series2)));

}

TEST_LIST = {
    { "testing_search", testing_search },
    { NULL, NULL }
};
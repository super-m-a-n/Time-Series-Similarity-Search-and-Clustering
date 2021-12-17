// file : assist_functions.cpp
#include "assist_functions.hpp"

//computes floor(log_2(n)) + 1
unsigned int get_lg(unsigned int n){
    int r = 0;

    while ( n > 0){
        n = n >> 1;
        r += 1;
    }

    return r;
}

//returns largest power of 2 smaller than given number
unsigned int largest_power_of_2_smaller_than(unsigned int n)
{
	// returns the number with the leftmost set bit of number set, and the rest all unset
	n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    return n ^ (n >> 1);
}
//Functions used for assistance with minor functionalities of the program
#pragma once
#include <vector>

//computes floor(log_2(n)) + 1
unsigned int get_lg(unsigned int n);

//Used for the filtering of the input curve (Only used for continuous frechet)
std::vector<float> filter_input_curve(std::vector<float> input_data);


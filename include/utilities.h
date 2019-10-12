// remove this one
#pragma once
#include <vector>
#include <cmath>
#include <iostream>

double log_base_(double num, double base);
void print_array (double * array, int n, int m, std::string msg);

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v);

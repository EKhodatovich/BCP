#pragma once
#include <vector>


double simpson(const std::vector<double>& f_values, const std::vector<double>& f_half_values,  double dx);


std::vector<double> evalFunc(double a, double b, unsigned int N, double func(std::vector<double>), double x, unsigned int m);
double bessFunc(double t, std::vector<double> params);

std::vector<double> BessIntegral(const std::vector<double> &grid, unsigned int m);
double BessIntegral(double x, std::vector<double> params);

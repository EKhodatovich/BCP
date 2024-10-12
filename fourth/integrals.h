#pragma once
#include <vector>


double trapezoid(const std::vector<double>& f_values,  double dx);

double simpson(const std::vector<double>& f_values, const std::vector<double>& f_half_values,  double dx);


std::vector<double> evalFunc(double a, double b, unsigned int N, double func(double, double, unsigned int), double x, unsigned int m);

double bessFunc(double t, double x, unsigned int m);
std::vector<double> evalBessWithSimpson(double x_a, double x_b, unsigned int m, unsigned int N, unsigned int N_t);
std::vector<double> evalBessWithTrapezoid(double x_a, double x_b, unsigned int m, unsigned int N, unsigned int N_t);

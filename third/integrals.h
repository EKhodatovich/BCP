#pragma once
#include <vector>


double func_(double x);

double leftRect(const std::vector<double>& f_values, double dx);

double rightRect(const std::vector<double>& f_values,  double dx);

double trapezoid(const std::vector<double>& f_values,  double dx);

double mean(const std::vector<double>& f_half_values,  double dx);

double simpson(const std::vector<double>& f_values, const std::vector<double>& f_half_values,  double dx);


std::vector<double> evalFunc(double a, double b, unsigned int N, double func(double));

double errFunc( double t);

double errInt(double x);

double bessFunc(double t, double x, double m);

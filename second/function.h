#pragma once


#include <vector>
class params
{
public:
    double m, a, U0, h;
};

double cot (double x);
double f_ (double x, const std::vector<double> &pars);
double df_ (double x, const std::vector<double> &pars);

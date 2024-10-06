#include "integrals.h"
#include <cmath>
#include <numbers>

double func_(double x)
{
    return 1 / (1 + x*x);
}


double leftRect(const std::vector<double>& f_values, const double dx)
{
    if (f_values.size() < 2) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    for (auto it = f_values.begin(); it != f_values.end() - 1; ++it)
    {
        result += (*it) * dx;
    }
    return result;
}

double rightRect(const std::vector<double>& f_values, const double dx)
{
    if (f_values.size() < 2) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    for (auto it = f_values.begin() + 1; it != f_values.end(); ++it)
    {
        result += (*it) * dx;
    }
    return result;
}

double trapezoid(const std::vector<double>& f_values, const double dx)
{
    if (f_values.size() < 2) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    for (auto it = f_values.begin(); it != f_values.end() - 1; ++it)
    {
        result += (*it + *(it+1)) / 2 * dx;
    }
    return result;
}

double mean(const std::vector<double>& f_half_values, const double dx)
{
    if (f_half_values.size() < 1) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    for (auto it = f_half_values.begin(); it != f_half_values.end(); ++it)
    {
        result += (*it) * dx;
    }
    return result;
}

double simpson(const std::vector<double>& f_values, const std::vector<double>& f_half_values, const double dx)
{
    if (f_half_values.size() < 1) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    auto half_it = f_half_values.begin();
    for (auto it = f_values.begin(); it != f_values.end()-1; ++it)
    {
        result += (*it + *(it + 1) + 4 * (*half_it)) / 6 * dx;
        ++half_it;
    }
    return result;
}


std::vector<double> evalFunc(double a, double b, unsigned int N, double func(double))
{
    std::vector<double> result;
    auto dx = (b-a)/N;
    for (unsigned int i = 0; i <=N; ++i)
    {
        double x = a + dx * i;
        result.push_back(func(x));
    }
    return result;
}

double errFunc(const double t)
{
    return std::exp(-(t*t)) * 2 / std::sqrt(std::numbers::pi);
}

double errInt(double x)
{
    int N = 10;
    double dx = 1.0/ N;
    return simpson(evalFunc(0, x, N, errFunc), evalFunc(dx, x-dx, N-1, errFunc), dx);
}

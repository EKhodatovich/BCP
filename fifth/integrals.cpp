#include "integrals.h"
#include <cmath>
#include <iostream>
#include <numbers>

using namespace std::numbers;


double simpson(const std::vector<double>& f_values, const std::vector<double>& f_half_values, const double dx)
{
    if (f_values.empty() || f_half_values.empty()) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    auto half_it = f_half_values.begin();
    for (auto it = f_values.begin(); it != f_values.end()-1; ++it)
    {
        result += (*it + *(it + 1) + (*half_it) * 4) / 6 * dx;
        ++half_it;
    }
    return result;
}

std::vector<double> evalFunc(double t_a, double t_b, unsigned int N, double func(double, std::vector<double>), std::vector<double> params)
{
    std::vector<double> result;
    auto dt = (t_b-t_a)/N;
    for (unsigned int i = 0; i <=N; ++i)
    {
        double t = t_a + dt * i;
        result.push_back(func(t, params));
    }
    return result;
}

double bessFunc(double t, std::vector<double>params)
{
    if (params.size() != 2) {
        std::cerr << "signature of bessFunc: (t, {m, x})";
    }
    unsigned int m = params[0];
    double x = params [1];
    return std::cos(m * t - x * std::sin(t))/pi;
}

std::vector<double> BessIntegral(const std::vector<double>& grid, unsigned int m)
{
    std::vector<double> result;
    for (auto x : grid)
    {
        result.push_back(BessIntegral(x, {static_cast<double>(m)}));
    }
    return result;
}

double BessIntegral(double x, std::vector<double> params)
{
    if (params.size() != 1) {
        std::cerr << "signature of BessIntegral: (x, {m})";
    }
    double m = params[0];
    const auto t_N = 100;
    auto f_values = evalFunc(0, pi, t_N, bessFunc, {m, x});
    auto f_half_values = evalFunc(0 + pi/t_N/2, pi-pi/t_N/2, t_N-1, bessFunc, {m, x});
    return simpson(f_values, f_half_values, pi/t_N);
}

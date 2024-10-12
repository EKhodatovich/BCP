#include "integrals.h"
#include <cmath>
#include <numbers>

using namespace std::numbers;

double trapezoid(const std::vector<double>& f_values, const double dx)
{
    if (f_values.size() < 2) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    for (auto it = f_values.begin(); it != f_values.end() - 1; ++it) {
        result += (*it + *(it+1)) / 2 * dx;
    }
    return result;
}

double simpson(const std::vector<double>& f_values, const std::vector<double>& f_half_values, const double dx)
{
    if (f_values.empty() || f_half_values.empty()) {
        throw std::bad_array_new_length();
    }
    double result = 0.0;
    auto half_it = f_half_values.begin();
    for (auto it = f_values.begin(); it != f_values.end()-1; ++it) {
        result += (*it + *(it + 1) + (*half_it) * 4) / 6 * dx;
        ++half_it;
    }
    return result;
}

std::vector<double> evalFunc(double t_a, double t_b, unsigned int N, double func(double, double, unsigned int), double x, unsigned int m)
{
    std::vector<double> result;
    auto dt = (t_b-t_a)/N;
    for (unsigned int i = 0; i <=N; ++i) {
        double t = t_a + dt * i;
        result.push_back(func(t, x, m));
    }
    return result;
}

double bessFunc(double t, double x, unsigned int m)
{
    return std::cos(m * t - x * std::sin(t))/pi;
}

std::vector<double> evalBessWithSimpson(double x_a, double x_b, unsigned int m, unsigned int N, unsigned int N_t)
{
    auto dx = (x_b-x_a)/ N;
    std::vector<double> result;
    for (auto x = x_a; x < x_b; x+=dx) {

        auto f_values = evalFunc(0, pi, N_t, bessFunc, x, m);
        auto f_half_values = evalFunc(0 + pi/N_t/2, pi-pi/N_t/2, N_t-1, bessFunc, x, m);
        result.push_back(simpson(f_values, f_half_values, pi/N_t));
    }
    return result;
}

std::vector<double> evalBessWithTrapezoid(double x_a, double x_b, unsigned int m, unsigned int N, unsigned int N_t)
{
    auto dx = (x_b-x_a)/ N;
    std::vector<double> result;
    for (auto x = x_a; x < x_b; x+=dx) {
        auto f_values = evalFunc(0, pi, N_t, bessFunc, x, m);
        result.push_back(trapezoid(f_values, pi/N_t));
    }
    return result;
}

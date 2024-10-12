#include "lagrangepolynom.h"
#include <cmath>
#include <numbers>
#include <stdexcept>

LagrangePolynom::LagrangePolynom(const double x_a, const double x_b, const unsigned int N)
{
    if (x_a >= x_b || N <= 1) {
        throw std::bad_array_new_length();
    }
    auto dx = (x_b - x_a ) / N;
    for (int i = 0; i < N; ++i) {
        grid_.push_back(0.5 * (x_a + x_b) + 0.5 * (x_b - x_a) * std::cos(std::numbers::pi * (2 * i + 1) / (2 * N)));
    }
}

double LagrangePolynom::evalFunc(double x)
{
    double result = 0.0;
    int i = 0;
    for (auto it : coefs_) {
        result+=(evalAuxPoly(x, i) * it);
        ++i;
    }
    return result;
}

void LagrangePolynom::fit(double func(double, std::vector<double>), const std::vector<double>& params)
{
    int i = 0;
    for (auto it : grid_) {
        coefs_.push_back(func(it, params) / evalAuxPoly(it, i));
        ++i;
    }
}

void LagrangePolynom::fit(const std::vector<double> &f_values)
{
    if (f_values.size() != grid_.size()) {
        throw std::invalid_argument("Vectors must be of the same size and contain at least two elements.");
    }
    int i = 0;
    for (auto it : grid_) {
        coefs_.push_back(f_values[i] / evalAuxPoly(it, i));
        ++i;
    }
}

double LagrangePolynom::evalAuxPoly(double x, int i)
{
    double result =1.0;
    int i_ = 0;
    if (i >= 0){
        for (auto it : grid_) {
            if (i!=i_) {
                result *= (x - it);
            }
            ++i_;
        }
    } else {
        for (auto it : grid_) {
            result *= (x -it);
        }
    }
    return result;
}

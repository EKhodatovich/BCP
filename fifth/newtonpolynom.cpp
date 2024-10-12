#include "newtonpolynom.h"
#include <cmath>
#include <stdexcept>

NewtonPolynom::NewtonPolynom(const double x_a, const double x_b, const unsigned int N)
{
    if (x_a >= x_b || N <= 1) {
        throw std::bad_array_new_length();
    }
    auto dx = (x_b - x_a ) / N;
    for (int i = 0; i < N; ++i) {
        grid_.push_back(0.5 * (x_a + x_b) + 0.5 * (x_b - x_a) * std::cos(std::numbers::pi * (2 * i + 1) / (2 * N)));
    }
}

void NewtonPolynom::fit(const std::vector<double> &f_values)
{
    if (grid_.size() != f_values.size() || grid_.size() < 2) {
        throw std::invalid_argument("Vectors must be of the same size and contain at least two elements.");
    }
    size_t N = grid_.size();
    coefs_.resize(N);
    std::vector<std::vector<double>> table(N, std::vector<double>(N));

    for (size_t i = 0; i < N; ++i) {
        table[i][0] = f_values[i];
    }

    for (size_t j = 1; j < N; ++j) {
        for (size_t i = 0; i < N - j; ++i) {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (grid_[i + j] - grid_[i]);
        }
    }

    for (size_t i = 0; i < N; ++i) {
        coefs_[i] = table[0][i];
    }
}

void NewtonPolynom::fit(double func(double, std::vector<double>), const std::vector<double>& params)
{
    size_t N = grid_.size();
    coefs_.resize(N);
    std::vector<std::vector<double>> table(N, std::vector<double>(N));

    for (size_t i = 0; i < N; ++i) {
        table[i][0] = func(grid_[i], params);
    }

    for (size_t j = 1; j < N; ++j) {
        for (size_t i = 0; i < N - j; ++i) {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (grid_[i + j] - grid_[i]);
        }
    }

    for (size_t i = 0; i < N; ++i) {
        coefs_[i] = table[0][i];
    }
}

double NewtonPolynom::evalFunc(const double x)
{
    double result = coefs_[0];
    double term = 1.0;
    for (int k = 1; k < grid_.size(); ++k) {
        term *= (x - grid_[k-1]);
        result += coefs_[k] * term;
    }
    return result;
}

double NewtonPolynom::evalAuxPoly(const double x, const int k)
{
    double result = 1.0;

    for (auto it = grid_.begin(); it != grid_.begin() + k; ++it) {
        result *= (x - *it);
    }
    return result;
}

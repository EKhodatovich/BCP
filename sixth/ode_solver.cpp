#include "ode_solver.h"
#include <iostream>
#include <limits>

ODE_Solver::ODE_Solver(double (*func)(double, double), double x0, double t0, double t_end, double step)
    : func_(func), x0_(x0), t0_(t0), t1_(t_end), step_(step) {}

std::vector<double> ODE_Solver::euler() {
    double t_ = t0_;
    double x_ = x0_;
    std::vector<double> result;
    while (t_ <= t1_) {
        result.push_back(x_);
        x_ += step_ * func_(t_, x_);
        t_ += step_;
    }
    return result;
}

std::vector<double> ODE_Solver::runge_kutta_2(const double alpha) {
    if (alpha < std::numeric_limits<double>::epsilon()) {
        std::cerr << "alpha is not positive" << std::endl;
    }

    double t_ = t0_;
    double x_ = x0_;
    std::vector<double> result;
    while (t_ <= t1_) {
        result.push_back(x_);
        double f = func_(t_, x_);
        double k1 = (1-alpha) * f;
        double k2 = alpha * func_(t_ + step_ / 2/alpha , x_ + step_ / 2/alpha * f);
        x_ += step_  * (k1 + k2);
        t_ += step_;
    }
    return result;
}

std::vector<double> ODE_Solver::runge_kutta_4() {
    double t_ = t0_;
    double x_ = x0_;
    std::vector<double> result;
    while (t_ <= t1_) {
        result.push_back(x_);
        double k1 = func_(t_, x_);
        double k2 = func_(t_ + step_ / 2, x_ + step_ / 2 * k1);
        double k3 = func_(t_ + step_ / 2, x_ + step_ / 2 * k2);
        double k4 = func_(t_ + step_, x_ + step_ * k3);
        x_ += (step_ / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        t_ += step_;
    }
    return result;
}

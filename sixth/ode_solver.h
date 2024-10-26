#pragma once

#include <vector>
class ODE_Solver {
public:
    ODE_Solver(double (*func)(double, double), double x0, double t0, double t_end, double step);

    std::vector<double> euler();
    std::vector<double> runge_kutta_2(double alpha);
    std::vector<double> runge_kutta_4();

private:
    double (*func_)(double, double);
    double x0_;
    double t0_;
    double t1_;
    double step_;
};


#pragma once

#include <vector>
class LagrangePolynom
{
public:
    LagrangePolynom(double x_a, double x_b, unsigned int N);
    double evalFunc(double x);
    void fit(double func(double, std::vector<double>), const std::vector<double>& params);
    void fit(const std::vector<double>& f_values);
    const std::vector<double>& getGrid() {return grid_;};
private:
    std::vector<double> grid_;
    std::vector<double> coefs_;

    double evalAuxPoly(double x, int i);

};

#include "function.h"
#include <cmath>
#include <exception>

double cot (const double x)
{
    return 1/ tan(x);
}

double f_ (const double x, const std::vector<double>& pars)
{
    if (pars.size()<4) {
        throw std::bad_exception();
    }
    double m = pars[0], a = pars[1], U0 = pars[2], h = pars[3];
    return sqrt(1/x - 1) - cot(  sqrt (2 * m * a*a * U0 / h/h * (1 - x))  );
}

double df_ (const double x, const std::vector<double> &pars)
{
    if (pars.size()<4) {
        throw std::bad_exception();
    }
    double m = pars[0], a = pars[1], U0 = pars[2], h = pars[3];
    double c = 2 * m * a*a * U0 / h/h ;
    return -1 / (2 * x*x * sqrt(1/x - x)) - \
            sqrt(c)/ (  2 * sqrt (1-x) * std::pow( sin( sqrt(c * (1-x)) ) , 2)    );
}

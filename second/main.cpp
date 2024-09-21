#include "function.h"
#include "mainwindow.h"

#include <QApplication>
#include <iostream>

#define pair std::pair<double,double>

#define EPS 1.0e-10

static const std::vector<double> pars {
    9.1e-28,
    1e-12,
    5e7 * 1.6e-12,
    1.054e-27
};

pair dichotomy_step(pair section)
{
    static std::unordered_map<double,double> f_values;
    if (f_values.empty())
    {
        f_values.try_emplace(section.first, f_(section.first, pars));
        f_values.try_emplace(section.second, f_(section.second, pars));
    }

    pair result;

    double avg = (section.first + section.second) / 2;
    double f_avg = f_(avg, pars);
    f_values.try_emplace(avg, f_avg);

    if (f_values.find(section.first)->second * f_avg < 0)
    {
        result.first = section.first;
        result.second = avg;
    }
    else
    {
        result.first = avg;
        result.second = section.second;
    }
    return result;
}

double dichotomy (double a, double b, double f(double, const std::vector<double>&))
{
    pair section = {a, b};
    unsigned int count = 0;
    while ( f((section.first + section.second) / 2, pars) !=0 && fabs(section.first - section.second) > EPS)
    {
        section = dichotomy_step(section);
        ++count;
    }
    std::cout << "f(" << (section.first + section.second) / 2 << ") = " << f((section.first + section.second) / 2, pars) << std::endl;
    std::cout << "count: " << count <<std::endl;
    return (section.first + section.second) / 2;
}


double simple_iterations (double x0,\
                         double f(double, const std::vector<double>&), \
                         double df(double, const std::vector<double>&)  )
{

    double f_value = f(x0, pars);
    //std::cout << "df = " << df(x0, pars) << std::endl;
    double lambda = 1. / df(x0, pars);

    double x_new = x0 - lambda * f_value;
    unsigned int count = 0;
    while (fabs (x_new - x0) > EPS)
    {
        x0 = x_new;
        f_value = f(x0, pars);
        x_new = x0 - lambda * f_value;
        ++count;
    }
    std::cout << "f(" << x_new << ") = " << f(x_new, pars) << std::endl;
    std::cout << "count: " << count <<std::endl;
    return x_new;
}

double newtons_method (double x0,\
                         double f(double, const std::vector<double>&), \
                         double df(double, const std::vector<double>&)  )
{

    double f_value = f(x0, pars);
    //std::cout << "df = " << df(x0, pars) << std::endl;
    double lambda = 1. / df(x0, pars);

    double x_new = x0 - lambda * f_value;
    unsigned int count = 0;
    while (fabs (x_new - x0) > EPS)
    {
        x0 = x_new;
        f_value = f(x0, pars);
        x_new = x0 - 1/df(x0, pars) * f_value;
        ++count;
    }
    std::cout << "f(" << x_new << ") = " << f(x_new, pars) << std::endl;
    std::cout << "count: " << count <<std::endl;
    return x_new;
}


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    //MainWindow w;
    //w.show();
    std::cout.precision(16);
    std::vector<double> params = {
        9.1e-28,
        1e-12,
        5e7 * 1.6e-12,
        1.054e-27};
    std::cout << "dichotomy: " <<std::endl << dichotomy(0.01, 0.25, &f_) << std::endl;
    std::cout << "simple iterations: " << std::endl << simple_iterations(0.13, &f_, &df_) << std::endl;
    std::cout << "newtons method: " << std::endl << newtons_method(0.13, &f_, &df_) << std::endl;

    return 0;//a.exec();
}


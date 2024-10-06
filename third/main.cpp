#include "mainwindow.h"
#include "integrals.h"
#include <QApplication>
#include <iostream>


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    MainWindow w;
    //w.show();
    const double a = -1.0;
    const double b = 1.0;
    const unsigned int N = 10;
    const auto f_values = evalFunc(a, b, N, func_);
    const auto dx = (b-a)/N;
    const auto f_half_values = evalFunc(a + dx/2 , b - dx/2, N-1, func_);

    const double anValue = std::numbers::pi/2.0;

    std::cout.precision(16);
    std::cout << "N: " << N << std::endl;
    std::cout << "Formula of left rectangles: " << leftRect(f_values, dx) - anValue << std::endl;
    std::cout << "Formula of right rectangles: " << rightRect(f_values, dx) - anValue << std::endl;
    std::cout << "Formula of trapezoid: " << trapezoid(f_values, dx) - anValue << std::endl;
    std::cout << "Formula of mean: " << mean(f_half_values, dx) - anValue << std::endl;
    std::cout << "Formula of Simpson: " << simpson(f_values, f_half_values, dx) - anValue << std::endl << std::endl;


    std::cout << "Err function: " << errInt(1.0) << std::endl;
    return app.exec();
}

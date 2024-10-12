#include "integrals.h"
#include "mainwindow.h"
#include <iostream>
#include <QApplication>
#include <QLineSeries>
#include <QtCharts>
#include <cmath>

//#include <boost/math/special_functions/bessel.hpp>
using namespace std::numbers;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
/*
    unsigned int N = 10000;
    auto dx = 2.0 * pi / N;
    auto besselM1 = evalBessWithSimpson(0, 2 * pi, 1, N);
    auto besselM0 = evalBessWithSimpson(0, 2 * pi, 0, N);


    std::vector<double>derBesselM0;
    for (auto it = besselM0.begin()+2; it != besselM0.end()-2; ++it) {
        derBesselM0.push_back(
            ( *(it-2) / 12
             -*(it-1) * 2/3
             +*(it+1) * 2/3
             -*(it+2) /12
             ) /dx);
    }
    std::vector <double> result;
    std::transform(
        besselM1.begin()+2,
        besselM1.end(),
        derBesselM0.begin(),
        std::back_inserter(result),
        [](double a, double b)
        {
            return (b+a);
        });


    besselM1 = evalBessWithTrapezoid(0, 2 * pi, 1, N);
    besselM0 = evalBessWithTrapezoid(0, 2 * pi, 0, N);


    derBesselM0.clear();
    for (auto it = besselM0.begin()+2; it != besselM0.end()-2; ++it) {
        derBesselM0.push_back(
            ( *(it-2) / 12
             -*(it-1) * 2/3
             +*(it+1) * 2/3
             -*(it+2) /12
             ) /dx);
    }

    std::vector <double> result_2;
    std::transform(
        besselM1.begin()+2,
        besselM1.end(),
        derBesselM0.begin(),
        std::back_inserter(result_2),
        [](double a, double b)
        {
            return (b+a);
        });

    QLineSeries series_1, series_2;
    for (int i=0; i<result.size()-2; ++i)
    {
        series_1.append({dx*i, result[i]});
        series_2.append({dx*i, result_2[i]});
    }

    QChart chart{};
    chart.addSeries(&series_1);
    chart.addSeries(&series_2);

    // const auto t_N = 10;
    // auto f_values = evalFunc(0, pi, t_N, bessFunc, pi, 0);
    // QLineSeries series;
    // for (int i=0; i<f_values.size()-2; ++i)
    // {
    //     series.append({0 + pi/10*i, f_values[i]});
    //     std::cout << f_values[i] << std::endl;
    // }
    // chart.addSeries(&series);


    chart.legend()->markers(&series_1).first()->setLabel("Simpson");
    chart.legend()->markers(&series_2).first()->setLabel("Trapezoid");

    QValueAxis axisX;
    QValueAxis axisY;
    axisY.setTickCount(21);
    //axisY.setRange(-3e-11, 3e-11);

    chart.addAxis(&axisX, Qt::AlignBottom);
    chart.addAxis(&axisY, Qt::AlignLeft);

    series_1.attachAxis(&axisX);
    series_2.attachAxis(&axisX);
    axisY.setLabelFormat("%.1e");
    series_1.attachAxis(&axisY);
    series_2.attachAxis(&axisY);

    // series.attachAxis(&axisX);
    // series.attachAxis(&axisY);



    QChartView view(&chart);
    view.showMaximized();
*/
    MainWindow window{};

    std::cout << trapezoid({1/pi, 0, -1/pi}, pi/2) << std::endl;
    std::cout << simpson(evalFunc(0, pi, 2, bessFunc, pi, 1), evalFunc(pi/4, pi * 3 / 4, 1, bessFunc, pi, 1), pi/2) << std::endl;
    //std::cout << boost::math::sph_bessel(1, 2*pi/10) << std::endl;
    window.show();
    return app.exec();
}

//#include "mainwindow.h"

#include "integrals.h"
#include "lagrangepolynom.h"
#include <QApplication>
#include <iostream>
#include <QtCharts>
#include <QLineSeries>

int main(int argc, char *argv[])
{
    std::vector <QColor> colors = {    QColor {0, 255, 0}, // GREEN
                                  QColor {0, 0, 255}, // BLUE
                                  QColor (255, 0, 0), // RED
                                  QColor (12, 21, 140), // VULCAN
                                  QColor (255, 0, 255), // PINK
                                  QColor (0, 0, 0), // BLACK
                                  QColor (255, 255, 0), // YELLOW
                                  QColor (0, 255, 255), // CYAN
                                  QColor (0, 102, 51), // FOREST GREEN
                                  QColor (0, 128, 255)}; // ORANGE

    QApplication a(argc, argv);

    double x_a = 0.0, x_b = 10;//2*std::numbers::pi;
    std::cout.precision(16);

    std::vector<QLineSeries*> series;
    for (int N = 10; N <= 80; N+=10) {
        LagrangePolynom lp(x_a, x_b, N);
        auto grid = lp.getGrid();
        auto dx = (x_b-x_a) / N;


        auto f_values = BessIntegral(grid, 0);
        lp.fit(f_values);
        QLineSeries* seria = new QLineSeries;

        int i = 0;
        for (auto x : grid) {
            x+=dx/2;
            auto simpleBess = BessIntegral(x, {0});
            //auto simpleBess = f_values[i];
            auto lpBess = lp.evalFunc(x);
            seria->append(x, std::fabs(lpBess - simpleBess) == 0 ? std::numeric_limits<double>::epsilon() : std::fabs(lpBess - simpleBess));
            std::cout << std::scientific << lpBess << " - "<< simpleBess << " = " << (std::fabs(lpBess - simpleBess) == 0 ? std::numeric_limits<double>::epsilon() : std::fabs(lpBess - simpleBess)) << std::endl  ;
            i++;
        }

        // lp.fit(BessIntegral, {0});
        // QLineSeries* seria = new QLineSeries;
        // std::cout << "N = " << N << std::endl;
        // for (auto x : grid) {
        //     x+=dx/2;
        //     auto simpleBess = BessIntegral(x, {0});
        //     auto lpBess = lp.evalFunc(x);
        //     seria->append(x, std::fabs(lpBess - simpleBess) == 0 ? std::numeric_limits<double>::epsilon() : std::fabs(lpBess - simpleBess));
        //     std::cout << std::scientific << lpBess << " - "<< simpleBess << " = " << (std::fabs(lpBess - simpleBess) == 0 ? std::numeric_limits<double>::epsilon() : std::fabs(lpBess - simpleBess)) << std::endl  ;
        // }


        std::cout << "-----------------------------------------------" << std::endl;
        QPen pen = seria->pen();
        pen.setWidth(3);
        pen.setBrush(colors[(N /10 - 1) % colors.size()]);
        seria->setPen(pen);
        series.push_back(seria);
    }
    QChart chart;
    for (auto it = series.begin(); it != series.end(); ++it) {
        chart.addSeries(*it);
        chart.legend()->markers(*it).first()->setLabel(QString::number((it - series.begin()+ 1) * 10));
    }

    QValueAxis axisX;
//    QValueAxis axisY;
    axisX.setRange(0, 10);
//    QLogValueAxis axisX;
    QLogValueAxis axisY;
//    axisX.setBase(10);
    axisY.setRange(1e-17, 1e4);
    axisY.setBase(10);
    axisY.setLabelFormat("%.0e");

    chart.addAxis(&axisX, Qt::AlignBottom);
    chart.addAxis(&axisY, Qt::AlignLeft);

    for (auto it = series.begin(); it != series.end(); ++it) {
        (*it)->attachAxis(&axisX);
        (*it)->attachAxis(&axisY);
    }



    QChartView view(&chart);
    view.showMaximized();

//    MainWindow w;
//    w.show();
    return a.exec();
}

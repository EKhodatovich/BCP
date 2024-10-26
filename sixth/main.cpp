#include "ode_solver.h"
#include <array>
#include <QtCharts>
#include <iostream>

using namespace std;

// Пример функции для задачи Коши
double exampleFunction(double t, double x)
{
    return -std::sin(t); // dx/dt = -x/2; решение: x = x0 * exp(-t/2)
}

double exampleSolution (double x0, double t)
{
    return std::cos(t);
}

int main(int argc, char** argv) {

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

    QApplication app{argc, argv};
    double x0 = 1;   // Начальное значение
    double t0 = 0;   // Начальное время
    double t_end = 3; // Конечное время

    std::array<std::vector<double>, 5> results;
    std::array<QLineSeries*, 5> series{};
    for (int i = 0; i < series.size(); ++i) {
        series[i] = new QLineSeries;
        QPen pen = series[i]->pen();
        pen.setWidth(3);
        pen.setBrush(colors[i]);
        series[i]->setPen(pen);
    }

    std::array<QString, 6> names = {"Euler", "Modified Euler", "Changed Euler", "RK-2 a=0.75", "RK-4", "An. solution"};

    for (int N = 10; N < 1000; N+=10)
    {
        double step = 1.0 / N;
        ODE_Solver solver(exampleFunction, x0, t0, t_end, step);

        results[0] = std::move(solver.euler());
        results[1] = std::move(solver.runge_kutta_2(0.5));
        results[2] = std::move(solver.runge_kutta_2(1));
        results[3] = std::move(solver.runge_kutta_2(0.75));
        results[4] = std::move(solver.runge_kutta_4());

        for (int i = 0; i < results.size(); ++i)
        {
            for (auto it = results[i].begin(); it != results[i].end(); it++) {
                double t = t0 + step * (it - results[i].begin());
                *it = std::fabs (*it - exampleSolution(x0, t) + std::numeric_limits<double>::epsilon());
            }
        }
        for (int i = 0; i < series.size(); ++i) {
            series[i]->append(QPointF{static_cast<double>(N), *std::max_element(results[i].begin(), results[i].end())});
        }
    }

    QChart chart;
    for (auto it = series.begin(); it != series.end(); ++it) {
        chart.addSeries(*it);
        chart.legend()->markers(*it).first()->setLabel(names[it - series.begin()]);
    }

    //QValueAxis axisX;
    //QValueAxis axisY;
    auto *axisX = new QLogValueAxis();
    auto *axisY = new QLogValueAxis();

    //axisX.setRange(0, 10);
    axisY->setRange(1e-15, 1e-1);
    axisX->setBase(10);
    axisY->setBase(10);

    axisY->setLabelFormat("%.0e");

    chart.addAxis(axisX, Qt::AlignBottom);
    chart.addAxis(axisY, Qt::AlignLeft);

    for (auto *it = series.begin(); it != series.end(); ++it) {
        (*it)->attachAxis(axisX);
        (*it)->attachAxis(axisY);
    }
    QChartView view(&chart);
    view.showMaximized();

    return QApplication::exec();
}

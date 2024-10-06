#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "integrals.h"
#include <QtCharts>
#include <QChartView>
#include <iostream>
#include <numbers>

//typedef double (*integral) (std::vector<double>, double);

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QChart *chart = new QChart();
    chart->setTitle("Function");
    std::vector<QString> intNames = {"leftRect", "rightRect", "trapezoid", "mean", "Simpson"};

    std::array<QLineSeries*, 5> series =
        {
        new QLineSeries(),
        new QLineSeries(),
        new QLineSeries(),
        new QLineSeries(),
        new QLineSeries(),
    };
    const double anValue = std::numbers::pi/2.0;
    for (unsigned int N = 10; N < 1000; N+=10)
    {
        auto dx = 2. / N;
        auto f_values = evalFunc(-1.0, 1.0, N, func_);
        auto f_half_values = evalFunc(-1.0+dx/2.0, 1.0-dx/2.0, N-1, func_);
        series[0]->append(QPointF{static_cast<qreal>(N), std::fabs(leftRect(f_values, dx) - anValue)});
        series[1]->append(QPointF{static_cast<qreal>(N), std::fabs(rightRect(f_values, dx) - anValue)});
        series[2]->append(QPointF{static_cast<qreal>(N), std::fabs(trapezoid(f_values, dx) - anValue)});
        series[3]->append(QPointF{static_cast<qreal>(N), std::fabs(mean(f_half_values, dx) - anValue)});
        series[4]->append(QPointF{static_cast<qreal>(N), std::fabs(simpson(f_values, f_half_values, dx) - anValue)+std::numeric_limits<double>::epsilon()});
    }
    for (auto it = 4; it < 5; ++it)
    {
        chart->addSeries(series[it]);
        chart->legend()->markers(series[it]).first()->setLabel(intNames[it]);
    }

    auto axisX = new QLogValueAxis;
    axisX->setBase(10);
    chart->addAxis(axisX, Qt::AlignBottom);

    auto axisY = new QLogValueAxis;
    axisY->setBase(10);
    chart->addAxis(axisY, Qt::AlignLeft);
    for (auto it = 4; it < 5; ++it)
    {
        series[it]->attachAxis(axisX);
        series[it]->attachAxis(axisY);
    }


    //chart->axes(Qt::Horizontal).first()->setRange(2, 100);
    QChartView* chartView = new QChartView(chart);
    chartView->show();
}


MainWindow::~MainWindow()
{
    delete ui;
}

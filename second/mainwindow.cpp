#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "function.h"
#include <QtCharts>
#include <QChartView>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    std::vector<double> pars {
        9.1e-28,
        1e-12,
        5e7 * 1.6e-12,
        1.054e-27
    };

    QChart *chart = new QChart();
    chart->setTitle("Function");
    QLineSeries* series = new QLineSeries(chart);
    for (double x = 0.068; x<0.125; x+=0.001)
    {
        QPointF point{x, f_(x, pars)};
        series->append(point);
    }
    chart->addSeries(series);
    chart->createDefaultAxes();
    chart->axes(Qt::Vertical).first()->setRange(-0.5, 0.5);
    chart->axes(Qt::Horizontal).first()->setRange(0.111, 0.125);
    QChartView* chartView = new QChartView(chart);
    ui->gridLayout->addWidget(chartView, 1, 2);
}

MainWindow::~MainWindow()
{
    delete ui;
}

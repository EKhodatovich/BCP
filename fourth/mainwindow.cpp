#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "integrals.h"
#include <iostream>

using namespace std::numbers;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow), dx_(2.0 * pi / N_),
    series_1( new QLineSeries()),
    series_2( new QLineSeries()),
    chart_ (new QChart()),
    chartView_ ( new QChartView(chart_)),
    axisX( new QValueAxis),
    axisY (new QValueAxis)
    //chartView_ (new QChartView())
{
    ui->setupUi(this);

    axisY->setTickCount(21);
    //axisY.setRange(-3e-11, 3e-11);
    std::cout.precision(16);
    evalIntegrals();
    makeSeries();
    chartView_->showMaximized();
}

void MainWindow::evalIntegrals()
{
    std::cout << "================================================\n"
              << "N = " << N_ << ", N_t = " << N_t_ << std::endl;
    auto besselM1 = std::move(evalBessWithSimpson(0, 2 * pi, 1, N_, N_t_));
    auto besselM0 = std::move(evalBessWithSimpson(0, 2 * pi, 0, N_, N_t_));



    std::vector<double>derBesselM0;
    for (auto it = besselM0.begin()+2; it != besselM0.end()-2; ++it) {
        derBesselM0.push_back(
            ( *(it-2) / 12
             -*(it-1) * 2/3
             +*(it+1) * 2/3
             -*(it+2) /12
             ) /dx_);
        std::cout << *it << " - " << _j0(dx_ * (it - besselM0.begin())) << " = " << *it - _j0(dx_ * (it - besselM0.begin())) << std::endl;
    }
    result.clear();
    std::transform(
        besselM1.begin()+2,
        besselM1.end()-2,
        derBesselM0.begin(),
        std::back_inserter(result),
        [](double a, double b)
        {
            //std::cout << a << " + " << b << " = " << a+b << std::endl;
            return (b+a);
        });

    std::cout << "--------------------------------------------" << std::endl;

    besselM1 = std::move(evalBessWithTrapezoid(0, 2 * pi, 1, N_, N_t_));
    besselM0 = std::move(evalBessWithTrapezoid(0, 2 * pi, 0, N_, N_t_));

    derBesselM0.clear();
    for (auto it = besselM0.begin()+2; it != besselM0.end()-2; ++it)
    {
        derBesselM0.push_back(
            ( *(it-2) / 12
             -*(it-1) * 2/3
             +*(it+1) * 2/3
             -*(it+2) /12
             ) /dx_);
        std::cout << *it << " - " << _j0(dx_ * (it - besselM0.begin())) << " = " << *it - _j0(dx_ * (it - besselM0.begin())) << std::endl;
    }
    result_2.clear();
    std::transform(
        besselM1.begin()+2,
        besselM1.end()-2,
        derBesselM0.begin(),
        std::back_inserter(result_2),
        [](double a, double b)
        {
            //std::cout << a << " + " << b << " = " << a+b << std::endl;
            return (b+a);
        });
}

void MainWindow::makeSeries()
{
    series_1 = new QLineSeries();
    series_2 = new QLineSeries();
    axisX = new QValueAxis;
    axisY = new QValueAxis;
    //chart_ = new QChart();

    for (int i=0; i<result.size()-2; ++i)
    {
        series_1->append({dx_*i, result[i]});
        series_2->append({dx_*i, result_2[i]});
    }
    //chart_ = new QChart;
    chart_->addSeries(series_1);
    chart_->addSeries(series_2);

    // const auto t_N = 10;
    // auto f_values = evalFunc(0, pi, t_N, bessFunc, pi, 0);
    // QLineSeries series;
    // for (int i=0; i<f_values.size()-2; ++i)
    // {
    //     series.append({0 + pi/10*i, f_values[i]});
    //     std::cout << f_values[i] << std::endl;
    // }
    // chart.addSeries(&series);


    chart_->addAxis(axisX, Qt::AlignBottom);
    chart_->addAxis(axisY, Qt::AlignLeft);
    axisY->setLabelFormat("%.1e");
    axisY->setRange(
        std::min(
            *std::min_element(result.begin(), result.end()-2),
            *std::min_element(result_2.begin(), result_2.end()-2)),
        std::max(
            *std::max_element(result.begin(), result.end()-2),
            *std::max_element(result_2.begin(), result_2.end()-2)));

    series_1->attachAxis(axisX);
    series_2->attachAxis(axisX);
    series_1->attachAxis(axisY);
    series_2->attachAxis(axisY);

    chart_->legend()->markers(series_1).first()->setLabel("Simpson");
    chart_->legend()->markers(series_2).first()->setLabel("Trapezoid");

    //series.attachAxis(&axisX);
    //series.attachAxis(&axisY);
    //chartView_ = new QChartView(chart_);
    //chartView_->showMaximized();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_comboBox_currentTextChanged(const QString &arg1)
{
    N_ = arg1.toInt();
    dx_ = 2.0 * pi / N_;
    clear();
    evalIntegrals();
    makeSeries();
}

void MainWindow::clear()
{
    //chartView_->close();
    //delete chart_;
    chart_->removeAllSeries();
    chart_->removeAxis(axisX);
    chart_->removeAxis(axisY);
    delete axisX;
    delete axisY;
}


void MainWindow::on_comboBox_2_currentTextChanged(const QString &arg1)
{
    N_t_ = arg1.toInt();
    clear();
    evalIntegrals();
    makeSeries();
}


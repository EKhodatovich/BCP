#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "predatorprey.h"
#include <iostream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow),
    chart_ (new QChart()),
    chartView_ ( new QChartView(chart_)),
    axisX( new QValueAxis),
    axisY (new QValueAxis)
{
    ui->setupUi(this);

    a_ = static_cast<double>(ui->spinBox->value());
    b_ = static_cast<double>(ui->spinBox_2->value());
    c_ = static_cast<double>(ui->spinBox_3->value());
    d_ = static_cast<double>(ui->spinBox_4->value());
    step_ = static_cast<double>(ui->comboBox->currentText().toDouble());
    preds_ = static_cast<double>(ui->spinBox_5->value());
    preys_ = static_cast<double>(ui->spinBox_6->value());


    paintSeries(evalSeries());
}

MainWindow::~MainWindow()
{
    delete ui;
    clear();
    delete chart_;
    delete chartView_;
}

void MainWindow::on_spinBox_valueChanged(int arg1)
{
    a_ = static_cast<double>(arg1);
    reEval();
}

void MainWindow::on_spinBox_2_valueChanged(int arg1)
{
    b_ = static_cast<double>(arg1);
    reEval();
}

void MainWindow::on_spinBox_3_valueChanged(int arg1)
{
    c_ = static_cast<double>(arg1);
    reEval();
}

void MainWindow::on_spinBox_4_valueChanged(int arg1)
{
    d_ = static_cast<double>(arg1);
    reEval();
}

void MainWindow::on_comboBox_currentTextChanged(const QString &arg1)
{
    step_ = arg1.toDouble();
    reEval();
}

void MainWindow::reEval()
{
    clear();
    paintSeries(evalSeries());
}

void MainWindow::clear()
{
    // Not only removes, but deletes all series
    chart_->removeAllSeries();
    chart_->removeAxis(axisX);
    chart_->removeAxis(axisY);
    delete axisX;
    delete axisY;
}

void MainWindow::paintSeries(const std::vector<QLineSeries*>& series)
{
    if (!names_.empty() && series.size() != names_.size() || series.empty()) {
        throw std::invalid_argument("Bad series number: it must be > 0 and == number of names");
    }
    for (int i = 0; i < series.size(); ++i) {
        chart_->addSeries(series[i]);
        chart_->legend()->markers(series[i]).first()->setLabel(names_[i]);
    }

    axisX = new QValueAxis;
    axisY = new QValueAxis;
    //auto *axisX = new QLogValueAxis();
    //auto *axisY = new QLogValueAxis();

    //axisX.setRange(0, 10);
    //axisY->setRange(1e-15, 1e-1);
    //axisX->setBase(10);
    //axisY->setBase(10);

    axisY->setLabelFormat("%.0e");

    chart_->addAxis(axisX, Qt::AlignBottom);
    chart_->addAxis(axisY, Qt::AlignLeft);

    for (auto it : series) {
        it->attachAxis(axisX);
        it->attachAxis(axisY);
    }
    chartView_->showMaximized();
}

std::vector<QLineSeries *> MainWindow::evalSeries()
{
    PredatorPrey model(a_, b_, c_, d_);
    double t = 0.0;  // Начальное время
    double endTime = 20; // Конечное время
    double x = preys_; // Начальное количество жертв
    double y = preds_;  // Начальное количество хищников

    std::vector <double> preyPopulation, predatorPopulation, timePoints;

    std::vector<QLineSeries *> series;
    for (int i =0 ; i < (solutionFlag_ ? 2 : 1); ++i) {
        series.push_back(new QLineSeries);
    }


    while (t <= endTime) {
        preyPopulation.push_back(x);
        predatorPopulation.push_back(y);
        timePoints.push_back(t);

        if (solutionFlag_){
            series[0]->append(QPointF{t, x});
            series[1]->append(QPointF{t, y});
        } else {
            series[0]->append(QPointF{x,y});
        }
        rungeKutta(model, x, y, t, step_);
        t += step_;
    }

    // for (size_t i = 0; i < timePoints.size(); ++i) {
    //     std::cout << "Time: " << timePoints[i] << ", Prey: " << preyPopulation[i] << ", Predators: " << predatorPopulation[i] << std::endl;
    // }
    return series;
}

void MainWindow::on_spinBox_5_valueChanged(int arg1)
{
    preds_ = static_cast<double>(arg1);
    reEval();
}

void MainWindow::on_spinBox_6_valueChanged(int arg1)
{
    preys_ = static_cast<double>(arg1);
    reEval();
}


void MainWindow::on_checkBox_clicked()
{
    solutionFlag_ = true;
    names_.clear();
    names_.push_back("Predators");
    names_.push_back("Preys");
    reEval();
}


void MainWindow::on_checkBox_2_clicked()
{
    solutionFlag_ = false;
    names_.clear();
    names_.push_back("Phase trajectory");
    reEval();
}


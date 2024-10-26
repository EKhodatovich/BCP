#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtCharts>
#include <QMainWindow>
#include <boost/core/noncopyable.hpp>

namespace Ui {
class MainWindow;
}

class MainWindow final: public QMainWindow, private boost::noncopyable
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() final;

private slots:
    void on_spinBox_valueChanged(int arg1);

    void on_spinBox_2_valueChanged(int arg1);

    void on_spinBox_3_valueChanged(int arg1);

    void on_spinBox_4_valueChanged(int arg1);

    void on_comboBox_currentTextChanged(const QString &arg1);

    void on_spinBox_5_valueChanged(int arg1);

    void on_spinBox_6_valueChanged(int arg1);

    void on_checkBox_clicked();

    void on_checkBox_2_clicked();

private:
    Ui::MainWindow *ui;
    std::vector<QLineSeries*> evalSeries();
    void paintSeries(const std::vector<QLineSeries*>& series);
    void clear();
    void reEval();

    bool solutionFlag_ = true;
    std::vector<QString> names_ = {"Predators", "Preys"};

    QChart* chart_;
    QChartView* chartView_;
    QValueAxis* axisX;
    QValueAxis* axisY;

    double a_, b_, c_, d_, step_, preds_, preys_;
};

#endif // MAINWINDOW_H

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLineSeries>
#include <QtCharts>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_comboBox_currentTextChanged(const QString &arg1);

    void on_comboBox_2_currentTextChanged(const QString &arg1);

private:
    Ui::MainWindow *ui;
    void evalIntegrals();
    void makeSeries();
    void clear();

    std::vector <double> result, result_2;
    unsigned int N_ = 1000;
    unsigned int N_t_ = 100;

    double dx_;
    QChart* chart_;
    QChartView* chartView_;


    QLineSeries* series_1;
    QLineSeries* series_2;
    QValueAxis* axisX;
    QValueAxis* axisY;
};

#endif // MAINWINDOW_H

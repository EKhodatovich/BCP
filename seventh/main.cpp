#include <QCoreApplication>
#include <QtCharts>

#include "mainwindow.h"

using namespace std;



int main(int argc, char** argv) {

    // double alpha = 10; // Параметры модели
    // double beta = 2;
    // double delta = 2;
    // double gamma = 10;

    // PredatorPrey model(alpha, beta, delta, gamma);

    // double x = 40.0; // Начальное количество жертв
    // double y = 9.0;  // Начальное количество хищников
    // double t = 0.0;  // Начальное время
    // double dt = 0.01; // Шаг времени
    // double endTime = 10; // Конечное время

    // vector<double> preyPopulation;
    // vector<double> predatorPopulation;
    // vector<double> timePoints;

    // std::array<QLineSeries*, 2> series = {new QLineSeries, new QLineSeries};

    // while (t <= endTime) {
    //     preyPopulation.push_back(x);
    //     predatorPopulation.push_back(y);
    //     timePoints.push_back(t);

    //     series[0]->append(QPointF{t, x});
    //     series[1]->append(QPointF{t, y});

    //     rungeKutta(model, x, y, t, dt);
    //     t += dt;
    // }

    // // Вывод результатов
    // for (size_t i = 0; i < timePoints.size(); ++i) {
    //     cout << "Time: " << timePoints[i] << ", Prey: " << preyPopulation[i] << ", Predators: " << predatorPopulation[i] << endl;
    // }


    // std::array<QString, 2> names = {"Preys", "Predators"};
    // QChart chart;
    // for (auto it = series.begin(); it != series.end(); ++it) {
    //     chart.addSeries(*it);
    //     chart.legend()->markers(*it).first()->setLabel(names[it - series.begin()]);
    // }

    // auto *axisX = new QValueAxis;
    // auto *axisY = new QValueAxis;
    // //auto *axisX = new QLogValueAxis();
    // //auto *axisY = new QLogValueAxis();

    // //axisX.setRange(0, 10);
    // //axisY->setRange(1e-15, 1e-1);
    // //axisX->setBase(10);
    // //axisY->setBase(10);

    // //axisY->setLabelFormat("%.0e");

    // chart.addAxis(axisX, Qt::AlignBottom);
    // chart.addAxis(axisY, Qt::AlignLeft);

    // for (auto *it = series.begin(); it != series.end(); ++it) {
    //     (*it)->attachAxis(axisX);
    //     (*it)->attachAxis(axisY);
    // }
    // QChartView view(&chart);
    // view.showMaximized();
    QApplication app{argc, argv};
    MainWindow window;
    window.show();
    return QApplication::exec();
}

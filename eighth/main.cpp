#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <QtCharts>
#include <QApplication>
#include <QChartView>
#include <fstream>

#include <boost/qvm.hpp>


using namespace std;

struct UVec {
    double u_;
    double v_;
};

bool UVect_compare_by_u (const UVec& a, const UVec& b) {
    return a.u_ < b.u_;
}


bool UVect_compare_by_v (const UVec& a, const UVec& b) {
    return a.v_ < b.v_;
}

UVec du(UVec vec) {
    UVec result {0,0};
    result.u_ = 998 * vec.u_ + 1998 * vec.v_;
    result.v_ = -999 * vec.u_ - 1999 * vec.v_;
    return result;
}

UVec f_an (double x) {
    UVec result {0,0};
    result.u_ = 2 * exp(-x) + exp(-1000 * x);
    result.v_ = - exp(-x) - exp(-1000 * x);
    return result;
}

vector<UVec> explicit_euler(UVec u0Vec, double h, double N)
{
    vector<UVec> result;
    UVec uVec = u0Vec;
    for (int i = 0; i < N; ++i){
        result.push_back(uVec);
        auto uPrimeVec = du(uVec);
        uVec.u_ += h * uPrimeVec.u_;
        uVec.v_ += h * uPrimeVec.v_;
    }
    return result;
}

/*
for (int i = 0; i < N; ++i){
    result.push_back(uVec);
    auto duVec = du(uVec);
    UVec uMidVec {uVec.u_ + h/2 * duVec.u_, uVec.v_ + h/2 * duVec.v_ };
    auto uPrimeVec = du(uMidVec);
    uVec.u_ += h * uPrimeVec.u_;
    uVec.v_ += h * uPrimeVec.v_;
}
*/

vector<UVec> semi_implicit_euler(UVec u0Vec, double h, double N) {
    vector<UVec> result;
    UVec uVec = u0Vec;

    boost::qvm::mat< double, 2, 2 > A;
    A.a[0][0] = 998;
    A.a[0][1] = 1998;
    A.a[1][0] = -999;
    A.a[1][1] = -1999;
    boost::qvm::mat<double, 2, 2> E;
    boost::qvm::set_identity(E);
    auto A_inv = boost::qvm::inverse(E - h * A);

    for (int i = 0; i < N; ++i){
        result.push_back(uVec);
        auto duVec = du(uVec);
        UVec uMidVec {h/2 * duVec.u_, h/2 * duVec.v_ };
        boost::qvm::vec<double, 2> U = {uMidVec.u_, uMidVec.v_};
        auto uPrimeVec = du(uMidVec);
        uVec.u_ += (A_inv * U + U).a[0];
        uVec.v_ += (A_inv * U + U).a[1];
    }
    return result;
}

/*

implicit:

for (int i = 0; i < N; ++i){
    result.push_back(uVec);
    auto duVec = du(uVec);
    UVec uNextVec {uVec.u_ + h * duVec.u_, uVec.v_ + h * duVec.v_ };
    auto uPrimeVec = du(uNextVec);
    uVec.u_ += h * uPrimeVec.u_;
    uVec.v_ += h * uPrimeVec.v_;
}

double a11 = 1 - h * 998;
double a12 = -h * 1998;
double a21 = h * 999;
double a22 = 1 + h * 1999;

// Решение системы уравнений
// Используем методом Крамера или подстановкой
double denom = a11 * a22 - a12 * a21; // определитель

uVec.u_ += (uNextVec.u_ * a22 - uNextVec.v_ * a12) / denom;
uVec.v_ += (a11 * uNextVec.v_ - a21 * uNextVec.u_) / denom;

Вывод значений на каждом шаге
cout << n+1 << "\t" << x1_new << "\t" << x2_new << endl;

*/

vector<UVec> implicit_euler(UVec u0Vec, double h, double N) {
    vector<UVec> result;
    UVec uVec = u0Vec;

    boost::qvm::mat< double, 2, 2 > A;
    A.a[0][0] = 998;
    A.a[0][1] = 1998;
    A.a[1][0] = -999;
    A.a[1][1] = -1999;
    boost::qvm::mat<double, 2, 2> E;
    boost::qvm::set_identity(E);
    auto A_inv = boost::qvm::inverse(E - h * A);

    for (int i = 0; i < N; ++i) {
        result.push_back(uVec);
        auto duVec = du(uVec);
        UVec uNextVec {h * duVec.u_, h * duVec.v_ };
        boost::qvm::vec<double, 2> U = {uNextVec.u_, uNextVec.v_};
        uVec.u_ += (A_inv * U).a[0];
        uVec.v_ += (A_inv * U).a[1];
    }

    return result;
}

int main (int argc, char* argv[])
{
    QApplication app(argc, argv);

    UVec u0Vec {3, -2};

    array<unique_ptr<QLineSeries>, 3> uSeries {make_unique<QLineSeries>(), make_unique<QLineSeries>(), make_unique<QLineSeries>()};
    array<unique_ptr<QLineSeries>, 3> vSeries {make_unique<QLineSeries>(), make_unique<QLineSeries>(), make_unique<QLineSeries>()};
    array<QString, 3> names {"Explicit", "Semi-implicit", "Implicit"};

    std::ofstream UExplicitFile("u_explicit.csv");
    std::ofstream USemiImplicitFile("u_semi_implicit.csv");
    std::ofstream UImplicitFile("u_implicit.csv");
    std::ofstream VExplicitFile("v_explicit.csv");
    std::ofstream VSemiImplicitFile("v_semi_implicit.csv");
    std::ofstream VImplicitFile("v_implicit.csv");

    for (int N = 10; N < 500000; N*=2) {
        cout << "Computing for " << N << " ..."<< endl;

        double h = 1.0 / N;
        auto explicit_ = explicit_euler(u0Vec, h, N);
        auto semiImplicit_ = semi_implicit_euler(u0Vec, h, N);
        auto implicit_ = implicit_euler(u0Vec, h, N);

        vector<UVec> f_values;
        for (int i = 0; i < explicit_.size(); ++i ) {
            f_values.push_back(f_an(i * h));
        }

        transform(explicit_.begin(),
                       explicit_.end(),
                       f_values.begin(),
                       explicit_.begin(),
                       [](UVec a, UVec b)
                        {
                            return UVec{fabs(a.u_ - b.u_),
                                       fabs(a.v_ - b.v_)};
                        });
        transform(semiImplicit_.begin(),
                       semiImplicit_.end(),
                       f_values.begin(),
                       semiImplicit_.begin(),
                       [](UVec a, UVec b)
                        {
                            return UVec{fabs(a.u_ - b.u_),
                                       fabs(a.v_ - b.v_)};
                        });
        transform(implicit_.begin(),
                       implicit_.end(),
                       f_values.begin(),
                       implicit_.begin(),
                       [](UVec a, UVec b)
                        {
                            //cout << a.u_ << " against " << b.u_ << endl ;
                            return UVec{fabs(a.u_ - b.u_) + std::numeric_limits<double>::epsilon(),
                                       fabs(a.v_ - b.v_) + std::numeric_limits<double>::epsilon()};
                        });

        uSeries[0]->append(N, max_element(explicit_.begin(), explicit_.end(), UVect_compare_by_u)->u_);
        uSeries[1]->append(N, max_element(semiImplicit_.begin(), semiImplicit_.end(), UVect_compare_by_u)->u_);
        uSeries[2]->append(N, max_element(implicit_.begin(), implicit_.end(), UVect_compare_by_u)->u_);
        vSeries[0]->append(N, max_element(explicit_.begin(), explicit_.end(), UVect_compare_by_v)->v_);
        vSeries[1]->append(N, max_element(semiImplicit_.begin(), semiImplicit_.end(), UVect_compare_by_v)->v_);
        vSeries[2]->append(N, max_element(implicit_.begin(), implicit_.end(), UVect_compare_by_v)->v_);


        UExplicitFile << N << "," << max_element(explicit_.begin(), explicit_.end(), UVect_compare_by_u)->u_ << endl;
        USemiImplicitFile << N << "," << max_element(semiImplicit_.begin(), semiImplicit_.end(), UVect_compare_by_u)->u_ << endl;
        UImplicitFile << N << "," << max_element(implicit_.begin(), implicit_.end(), UVect_compare_by_u)->u_ << endl;
        VExplicitFile << N << "," << max_element(explicit_.begin(), explicit_.end(), UVect_compare_by_v)->v_ << endl;
        VSemiImplicitFile << N << "," << max_element(semiImplicit_.begin(), semiImplicit_.end(), UVect_compare_by_v)->v_ << endl;
        VImplicitFile << N << "," << max_element(implicit_.begin(), implicit_.end(), UVect_compare_by_v)->v_ << endl << endl;
    }

    if (!names.empty() && uSeries.size() != names.size() || uSeries.empty()) {
        throw invalid_argument("Bad series number: it must be > 0 and == number of names");
    }
    auto uChart = new QChart();
    auto vChart = new QChart();
    for (int i = 0; i < uSeries.size(); ++i) {
        uChart->addSeries(uSeries.at(i).get());
        uChart->legend()->markers(uSeries.at(i).get()).first()->setLabel(names.at(i));
        vChart->addSeries(vSeries.at(i).get());
        vChart->legend()->markers(vSeries.at(i).get()).first()->setLabel(names.at(i));
    }

    // auto axisX = new QValueAxis;
    // auto axisY = new QValueAxis;
    auto *uAxisX = new QLogValueAxis;
    auto *uAxisY = new QLogValueAxis;

    auto *vAxisX = new QLogValueAxis;
    auto *vAxisY = new QLogValueAxis;

    // uAxisX->setRange(10000, 100000);
    // uAxisX->setRange(10000, 100000);
    uAxisY->setRange(1e-6, 1e+4);
    vAxisY->setRange(1e-6, 1e+4);

    uAxisX->setBase(10);
    uAxisY->setBase(10);
    vAxisX->setBase(10);
    vAxisY->setBase(10);

    uAxisY->setLabelFormat("%.0e");
    vAxisY->setLabelFormat("%.0e");

    uChart->addAxis(uAxisX, Qt::AlignBottom);
    uChart->addAxis(uAxisY, Qt::AlignLeft);
    vChart->addAxis(vAxisX, Qt::AlignBottom);
    vChart->addAxis(vAxisY, Qt::AlignLeft);

    for (auto it =uSeries.begin(); it != uSeries.end(); ++it) {
        (*it)->attachAxis(uAxisX);
        (*it)->attachAxis(uAxisY);
    }
    for (auto it =vSeries.begin(); it != vSeries.end(); ++it) {
        (*it)->attachAxis(vAxisX);
        (*it)->attachAxis(vAxisY);
    }

    // uChart->createDefaultAxes();
    auto* uChartView = new QChartView(uChart);
    uChartView->setWindowTitle("U");
    uChartView->showMaximized();

    auto* vChartView = new QChartView(vChart);
    vChartView->setWindowTitle("V");
    vChartView->showMaximized();

    return app.exec();
}

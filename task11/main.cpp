#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <fstream>

using namespace std;
using namespace boost::numeric;

static const double eigenValue = 1/2.0;

double Psi_exact (double x) {
    return exp( -x*x / 2.0 );
}

double U(double x) {
    return 1/2.0 * x*x;
}


void tridiagonalSolve(const vector<double>& a, const vector<double>& b, const vector<double>& c,
                      const ublas::vector<double>& d, ublas::vector<double>& y, int n) {
    vector<double> alpha(n), beta(n);
    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * alpha[i - 1];
        alpha[i] = c[i] / denom;
        beta[i] = (d[i] - a[i] * beta[i - 1]) / denom;
    }
    y[n - 1] = beta[n-1];

    for (int i = n - 2; i >= 0; --i) {
        y[i] = beta[i] - alpha[i] * y[i + 1];
    }
}


std::pair<double, ublas::vector<double>> SolveShred(const std::vector<double>& grid, double (*U)(double), double precision) {
    const int n = grid.size();
    const double h = grid[1] - grid[0];

    double main_diag_value = 1 / (h * h);
    double off_diag_value = -1 / (2 * h * h);

    vector<double> a(n, off_diag_value), b(n, main_diag_value), c(n, off_diag_value), d(n);
    for (int i = 0; i < n; ++i) {
        b[i] += U(grid[i]);
    }

    // Инициализация начального собственного вектора
    ublas::vector<double> eigenvector(n, 1.0);
    for (int i = 0; i < n; ++i) {
        eigenvector[i] += 1.0 * i/n;
    }

    // Начальное значение собственного числа
    double eigenValue = 0.0;
    double prevEigenValue = 0.0;
    double prevNorm = ublas::norm_2(eigenvector);
    int numberOfSteps = 0;

    while (true) {
        ublas::vector<double> nextEigenvector(n);
        tridiagonalSolve(a, b, c, eigenvector, nextEigenvector, n);

        auto norm = ublas::norm_2(nextEigenvector);
        eigenValue = norm / prevNorm;
        if (std::abs(eigenValue - prevEigenValue) < precision) {
            eigenvector = nextEigenvector;
            prevNorm = norm;
            break;
        }
        eigenvector = nextEigenvector;
        prevNorm = norm;
        prevEigenValue = eigenValue;
        numberOfSteps++;
    }
    cout << "number of steps: " << numberOfSteps << endl;
    eigenvector /= prevNorm;
    return {eigenValue, eigenvector};
}

int main()
{
    string strongDelim = "========================================\n";
    string delim =       "----------------------------------------\n";

    cout.precision(10);
    for (int precPower = 3; precPower <= 10; ++precPower)
    {
        double precision = pow(10, -precPower);
        cout << strongDelim << strongDelim << "precision = " << precision << endl << strongDelim << strongDelim;

        ofstream ff;
        vector<double> maxValues;
        vector<int> NValues;

        for (int N = 100; N <= 1000000; N *=2)
        {
            cout << delim << "N_x = " << N << endl << delim;
            double x1 = -10;
            double x2 = 10;
            double h = (x2-x1) / N;

            vector<double> grid;
            for (int i = 0; i < N; ++i) {
                grid.push_back(x1 + h*i);
            }


            auto result = SolveShred(grid, U, precision);
            auto result_vector = result.second;
            auto result_value = result.first;

            ublas::vector<double> exact_vector (N);
            for (int i = 0; i < N; ++i)
            {
                exact_vector[i] = Psi_exact(grid[i]);
            }
            auto norm = ublas::norm_2(exact_vector);
            exact_vector /= norm;

            for (int i = 0; i < N; ++i) {
                double a = result_vector[i];
                double b = exact_vector[i];
                // if (i % (N/ 100) == 0) {
                //     cout << "x = " << grid[i] << "  :\t" << a << ", " << b << endl;
                // }
                result_vector[i] = fabs(a - b);
            }
            double max = *max_element(result_vector.begin(), result_vector.end());
            cout << "max: " << max << endl;
            cout << "eigen value: " << 1.0/result_value << endl;
            maxValues.push_back(max);
            NValues.push_back(N);
        }

        std::string filename = "Prec" + std::to_string( precPower) + ".csv";
        ff.open(filename);
        for (int i = 0; i < maxValues.size(); ++i) {
            ff << NValues[i] << "," << maxValues[i] << endl;
        }
        ff.close();
    }
    return 0;
}

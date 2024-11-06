#include <iostream>
#include <cmath>
#include <numbers>
#include <vector>
#include <fstream>

using namespace std;

// y'' = f (x) = cos(x)
double f(double x) {
    return cos(x);
}

// General solution is
// y = -cos(x) + Ax + B

double y_an(const double x, const double A, const double B) {
    return -cos(x) + A * x + B;
}

void coeffForDirichDirich(const double y_left, const double y_right, double& A, double& B) {
    B = 1.0 / 2 * (y_right + y_left);
    A = 1.0 / numbers::pi * (y_right - y_left);
}

void coeffFor2BoundCondOnDeriv(const double y_left, const double y_right, double& A, double& B) {
    if (y_right - y_left - 2 > numeric_limits<double>::epsilon()) {
        throw invalid_argument("For two boundary conditions on the derivative of a function, the values ​​must differ strictly by 2");
    }
    A = y_left + 1;
    // B could be any - considered to be 0
    B = 0;
}

void coeffForNeumDirich(const double y_left, const double y_right, double& A, double& B) {
    A = y_left + 1;
    B = y_right - A*numbers::pi/2.0;
}

void coeffForDirichNeum(const double y_left, const double y_right, double& A, double& B) {
    A = y_right - 1;
    B = y_left + A*numbers::pi/2.0;
}

void tridiagonalSolve(const vector<double>& a, const vector<double>& b, const vector<double>& c,
                      const vector<double>& d, vector<double>& y, int n) {
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

// Решение для граничных условий на значения функции
vector<double> solveWithBoundaryValues(vector<double> grid, double A, double B, bool leftDirichle = true, bool rightDirichle = true) {
    if (!leftDirichle && !rightDirichle) {
        throw (invalid_argument("Both boundary conditions cannot be Neumanns'"));
    }

    int n = static_cast<int>(grid.size());

    vector<double> a(n, 1.0), b(n, -2.0), c(n, 1.0), d(n);

    if (leftDirichle) {
        b[0] = 1;
        c[0] = 0;
    } else {
        b[0] = -1;
        c[0] = 1;
    }
    if (rightDirichle) {
        b[n-1] = 1;
        a[n-1] = 0;
    } else {
        a[n-1] = -1;
        b[n-1] = 1;
    }

    double h = numbers::pi / n;
    for (int i = 0; i < n; ++i) {
        double x = grid[i];
        d[i] = h * h * cos(x);
    }

    if (leftDirichle) {
        d[0] = A;
    } else {
        d[0] = h * A;
    }

    if (rightDirichle) {
        d[n-1] = B;
    } else {
        d[n-1] = h * B;
    }

    vector<double> y(n);
    tridiagonalSolve(a, b, c, d, y, n);

    return (y);
}

int main() {

    ofstream ff("result_2.csv");

    string delim = "=======================================================\n";

    for (int n = 100; n < 10000000; n*=2)
    {
        cout << delim << "\t N = " << n << endl << delim;

        double h = numbers::pi / n;
        vector<double> grid(n);
        for (int i = 0; i < n; ++i) {
            grid[i] = -numbers::pi/2 + i * h;
        }

        // Пример 1: Граничные условия Дирихле-Дирихле
        double alpha = 1.0, beta = -1.0;
        auto values = solveWithBoundaryValues(grid, alpha, beta);

        double a = 0, b = 0;
        coeffForDirichDirich(alpha, beta, a, b);

        for (int i = 0; i < values.size(); ++i) {
            values[i] = fabs(values[i] - y_an(grid[i], a, b));
            if (i % (n/10) == 0) {
                cout << values [i] << endl;
            }
        }
        double max = *max_element(values.begin(), values.end());
        ff << n << "," << max;


        cout << delim;

        // Пример 2: Граничные условия Неймана-Дирихле
        alpha = 2;
        beta = 1;
        values = solveWithBoundaryValues(grid, alpha, beta, false, true);

        coeffForNeumDirich(alpha, beta, a, b);
        for (int i = 0; i < values.size(); ++i) {
            values[i] = fabs(values[i] - y_an(grid[i], a, b));
            if (i % (n/10) == 0) {
                cout << values [i] << endl;
            }
        }
        max = *max_element(values.begin(), values.end());
        ff << "," << max;


        cout << delim;

        // Пример 3: Граничные условия Дирихле-Неймана
        values = solveWithBoundaryValues(grid, alpha, beta, true, false);

        coeffForDirichNeum(alpha, beta, a, b);
        for (int i = 0; i < values.size(); ++i) {
            values[i] = fabs(values[i] - y_an(grid[i], a, b));
            if (i % (n/10) == 0) {
                cout << values [i] << endl;
            }
        }
        max = *max_element(values.begin(), values.end());
        ff << "," << max << endl;

    }

    return 0;
}

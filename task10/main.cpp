#include <cmath>
#include <numbers>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


/*
 * 2/l * integral ( 0; l; f(x) * sin (pi * n * x / l) dx)
 * for f(x) = sin (pi * x / l) will be
 * - 2 * sin (pi * n) / pi / (n**2 - 1)
 */

double FourierCoeff (int n, double l = 0) {
    if (n == 1) {
        return 1;
    }
    return - 2 * sin(numbers::pi * n) / numbers::pi / (n*n -1);
}


/*
 * U (x,t) = SUM ( n=0; inf) FourierCoeff(n) * sin (pi * n * x / l) * exp ( (pi * n / l)^2 * t)
 */

double u_an(double x, double t, double l) {
    double result = 0;
    for (int n = 1; n < 100; ++n) {

        double a = FourierCoeff(n);
        double b = sin (numbers::pi * n * x / l);
        double c = exp ( - pow((numbers::pi * n / l), 2) * t);
        result += a * b * c;
    }
    return result;
}

double timeBoundary(double x, double l) {
    return sin(numbers::pi * x / l);
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
vector<double> solveWithBoundaryValues(vector<double> prev_layer, double tau, double h) {

    int n = static_cast<int>(prev_layer.size());
    double r = tau /h/h / 2;

    vector<double>
        a(n, -r),
        b(n, 1 + 2*r),
        c(n, -r),
        d(n);
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = 0.0;

    a[n-1] = 0.0;
    b[n-1] = 1.0;
    c[n-1] = 0.0;

    /*
     * d__i = V_m_j + tau/2 * L[V_m_j]
     * L[V_m_j] = (V_m_j+1 - 2V_m_j + V_m_j-1) / h^2
     */

    for (int i = 1; i < d.size() - 1; ++i) {
        d[i] = prev_layer[i] + r * (prev_layer[i-1] - 2 * prev_layer[i] + prev_layer[i+1]);
    }
    d[0] = 0;
    d[n-1] = 0;

    vector<double> y(n);
    tridiagonalSolve(a, b, c, d, y, n);

    return (y);
}



int main()
{
    string strongDelim = "========================================\n";
    string delim =       "----------------------------------------\n";

    double maxTime = 1;
    for (int N_x = 100; N_x <= 10000; N_x *=10) {
        cout << strongDelim << strongDelim << "N_x = " << N_x << endl << strongDelim << strongDelim;

        ofstream ff;
        vector<double> tMaxValues;
        vector<int> NTValues;
        double h = 1.0 / N_x;
        for (int N_t = 100; N_t < 10000; N_t *=2) {
            cout << delim << "N_t = " << N_t << endl << delim;
            double tau = maxTime / N_t;

            vector<double> layer(N_x);
            for (int i = 0; i < N_x; ++i ) {
                double x = h * i;
                layer[i] = timeBoundary(x, 1);
            }

            for (int i = 0; i < N_t; ++i) {
                layer = solveWithBoundaryValues(layer, tau, h);
            }

            for (int i = 0; i < N_x; ++i) {
                double a = layer[i];
                double b = u_an( h*i, tau * (N_t -1), 1);
                if (i % (N_x/ 100) == 0) {
                    cout << a << ", " << b << endl;
                }
                layer[i] = fabs(a - b);
            }
            double max = *max_element(layer.begin(), layer.end());
            cout << max << endl;
            tMaxValues.push_back(max);
            NTValues.push_back(N_t);
        }

        std::string filename = "X_" + std::to_string(N_x) + ".csv";
        ff.open(filename);
        for (int i = 0; i < tMaxValues.size(); ++i) {
            ff << NTValues[i] << "," << tMaxValues[i] << endl;
        }
        ff.close();
    }
    return 0;
}

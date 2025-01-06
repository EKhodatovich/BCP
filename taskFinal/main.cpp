#include <algorithm>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>

constexpr double L = 10.0;       // Длина интервала
int N = 10000;           // Количество узлов
double T = 1.0/(N/100)/(N/100);        // Время моделирования
double dx = 2 * L / (N - 1); // Шаг по пространству
double dt = dx*dx / 100;      // Шаг по времени
const std::complex<double> i(0.0, 1.0); // Мнимая единица

// Параметры схемы Рунге-Кутты
const double lambda = 1;
const double c = 1;

using ComplexVec = std::vector<std::complex<double>>;

// Функция начальных условий
ComplexVec initialize(double lambda, double c, const std::vector<double>& xgrid ) {
    ComplexVec A(N);
    for (int j = 0; j < N; ++j) {
        double x = xgrid[j];
        A[j] = c * lambda / cosh(lambda * x);
    }
    return A;
}

// Вычисление второй производной с помощью центральных разностей
ComplexVec compute_second_derivative(const ComplexVec& A) {
    ComplexVec d2A(N);
    for (int j = 1; j < N - 1; ++j) {
        d2A[j] = (A[j + 1] - 2.0 * A[j] + A[j - 1]) / (dx * dx);
    }
    return d2A;
}

// Правая часть уравнения (дифференциальный оператор)
ComplexVec compute_rhs(const ComplexVec& A) {
    ComplexVec rhs(N);
    ComplexVec d2A = compute_second_derivative(A);
    for (int j = 1; j < N - 1; ++j) {
        rhs[j] = i * (2.0 * std::norm(A[j]) * A[j] + d2A[j]) * (-1.0);
    }
    return rhs;
}

// Метод Рунге-Кутты второго порядка
ComplexVec runge_kutta_step(const ComplexVec& A, double dt) {
    ComplexVec k1 = compute_rhs(A);
    ComplexVec A_mid(N);
    for (int j = 0; j < N; ++j) {
        A_mid[j] = A[j] + 0.5 * dt * k1[j];
    }
    ComplexVec k2 = compute_rhs(A_mid);
    ComplexVec A_next(N);
    for (int j = 0; j < N; ++j) {
        A_next[j] = A[j] + dt * k2[j];
    }
    return A_next;
}

// Запись модуля A(x, t) в файл
void save_to_file(const ComplexVec& A, double t, std::ofstream& file, const std::vector<double>& xgrid) {
    for (int j = 0; j < N; ++j) {
        double x = xgrid[j];
        file << x << " " << t << " " << std::abs(A[j]) << "\n";
    }
    file << "\n";
}

int main() {
    std::vector<double> xgrid(N);
    for (int j = 0; j < N; ++j)
    {
        xgrid.at(j) = -L + j * dx;
    }

    ComplexVec A = initialize(lambda, c, xgrid);
    std::ofstream file("errors.dat");
    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing.\n";
        return 1;
    }

    for (double t = 0.0; t < T; t += dt) {
        // save_to_file(A, t, file, xgrid);
        A = runge_kutta_step(A, dt);
    }
    ComplexVec an_solution = A;

    for (int n = 10; n < 10000; n*=2)
    {
        N = n;
        dx = 2 * L / (N - 1);
        dt = dx*dx / 100;

        std::vector<double> xgrid(N);
        for (int j = 0; j < N; ++j)
        {
            xgrid.at(j) = -L + j * dx;
        }

        ComplexVec A = initialize(lambda, c, xgrid);

        for (double t = 0.0; t < T; t += dt) {
            // save_to_file(A, t, file, xgrid);
            A = runge_kutta_step(A, dt);
        }
        std::vector<double> result(N);
        for (int j = 0; j < N; ++j) {
            result.at(j) = std::abs(A[j] - an_solution[j * 10000 / N]);
        }
        double max = *std::max_element(result.begin(), result.end());
        std::cout << N << " " << max << "\n";
        file << N << " " << max << std::endl;
    }
    std::cout << "\n";
    file.close();
    file = std::ofstream("time_errors.dat");
    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing.\n";
        return 1;
    }
    N = 10000;
    dx = 2 * L / (N - 1);
    dt = dx*dx / 100;
    for (int j = 0; j < N; ++j)
    {
        xgrid.at(j) = -L + j * dx;
    }

    for (int n = 100; n <= 1000; n+=10)
    {
        dt = dx*dx / n * 10;
        A = initialize(lambda, c, xgrid);
        for (double t = 0.0; t < T; t += dt) {
            // save_to_file(A, t, file, xgrid);
            A = runge_kutta_step(A, dt);
        }
        std::vector<double> result(N);
        for (int j = 0; j < N; ++j) {
            result.at(j) = std::abs(A[j] - an_solution[j]);
        }
        double max = *std::max_element(result.begin(), result.end());
        std::cout << 1000.0 / n << " " << max << std::endl;
        file << 1000.0 / n << " " << max << std::endl;
    }

    file.close();
    return 0;
}

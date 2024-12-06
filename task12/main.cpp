#include <iostream>
#include <numbers>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>

using namespace std;

// Параметры сигнала


// Функция для генерации сигнала
void generateSignal(vector<complex<double>>& signal, double a0, double a1, double w0, double w1, double T, int N) {
    for (int i = 0; i < N; ++i) {
        double t = T * i / (N - 1);
        double value = a0 * sin(w0 * t) + a1 * sin(w1 * t);
        signal[i] = complex<double>(value, 0.0);
    }
}

// Функция для применения окна Ханна
void applyHannWindow(vector<complex<double>>& signal, int N) {
    for (int i = 0; i < N; ++i) {
        double hann_value = 0.5 * (1 - cos(2 * numbers::pi * i / (N - 1)));
        signal[i] *= hann_value;
    }
}

// Функция для выполнения FFT
void FFT(vector<complex<double>>& x) {
    int N = x.size();
    if (N <= 1) return;

    vector<complex<double>> even(N / 2);
    vector<complex<double>> odd(N / 2);

    for (int i = 0; i < N / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    FFT(even);
    FFT(odd);

    for (int k = 0; k < N / 2; ++k) {
        complex<double> t = polar(1.0, -2 * numbers::pi * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

// Функция для вычисления спектра мощности
vector<double> powerSpectrum(const vector<complex<double>>& fft_result, int N) {
    vector<double> power(N / 2);
    for (int i = 0; i < N / 2; ++i) {
        power[i] = norm(fft_result[i]);
    }
    return power;
}

// Основная функция
int main() {
    const double a0 = 1.0;
    const double a1 = 1.0;
    const double w0 = 5.1;
    const double w1 = 25.1;
    const double T = 2 * numbers::pi;
    for (int power = 6; power < 20 ; power+=2) // количество точек дискретизации
    {
        ofstream ff_rect, ff_hann;
        ff_rect.open("N_" + to_string(power) + "_rect.csv");
        ff_hann.open("N_" + to_string(power) + "_hann.csv");
        int N = static_cast<int>(pow(2, power));

        const double sampling_rate = N / T;

        // Генерация сигнала
        vector<complex<double>> signal(N);
        generateSignal(signal, a0, a1, w0, w1, T, N);

        // Применение прямоугольного окна (без изменений сигнала)
        vector<complex<double>> signal_rect = signal;

        // Применение окна Ханна
        vector<complex<double>> signal_hann = signal;
        applyHannWindow(signal_hann, N);

        // Вычисление FFT для обоих сигналов
        FFT(signal_rect);
        FFT(signal_hann);

        // Вычисление спектра мощности
        vector<double> power_rect = powerSpectrum(signal_rect, N);
        vector<double> power_hann = powerSpectrum(signal_hann, N);

        cout << "Frequency\tRectangular Window\tHann Window\n";
        for (int i = 0; i < power_rect.size(); ++i) {
            double freq = i * sampling_rate / N ;
            if (i < 40) {
                cout << freq * 2 * numbers::pi << "\t" << power_rect[i] << "\t" << power_hann[i] << endl;
            }
            ff_rect << freq * 2 * numbers::pi << "," << power_rect.at(i) << endl;
            ff_hann << freq * 2 * numbers::pi << "," << power_hann.at(i) << endl;
        }

        ff_rect.close();
        ff_hann.close();
    }

    return 0;
}

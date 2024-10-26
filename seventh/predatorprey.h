#pragma once


struct PredatorPrey {
    double alpha; // Уровень размножения жертвы
    double beta;  // Уровень поглощения жертвой хищника
    double delta; // Уровень размножения хищника
    double gamma; // Уровень смертности хищника

    PredatorPrey(double a, double b, double d, double g)
        : alpha(a), beta(b), delta(d), gamma(g) {}

    void derivatives(double x, double y, double& dxdt, double& dydt) const {
        dxdt = alpha * x - beta * x * y; // Уравнение для жертвы
        dydt = delta * x * y - gamma * y; // Уравнение для хищника
    }
};

void rungeKutta(PredatorPrey& model, double& x, double& y, double t, double dt) {
    double k1x = 0, k1y = 0, k2x = 0, k2y = 0;

    model.derivatives(x, y, k1x, k1y);
    model.derivatives(x + k1x * dt / 2, y + k1y * dt / 2, k2x, k2y);

    x += dt * k2x;
    y += dt * k2y;
}

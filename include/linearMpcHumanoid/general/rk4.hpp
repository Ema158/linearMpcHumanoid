#pragma once

#include <Eigen/Dense>

template <typename DynamicsFunction>
Eigen::VectorXd rk4Step(
    DynamicsFunction f,
    const Eigen::VectorXd& x,
    double t,
    double dt)
{
    Eigen::VectorXd k1 = f(x, t);
    Eigen::VectorXd k2 = f(x + 0.5 * dt * k1, t + 0.5 * dt);
    Eigen::VectorXd k3 = f(x + 0.5 * dt * k2, t + 0.5 * dt);
    Eigen::VectorXd k4 = f(x + dt * k3, t + dt);

    return x + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}
#pragma once

#include <Eigen/Dense>

template <typename DynamicsFunction>
Eigen::VectorXd rk4Step(
    DynamicsFunction f,
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& u,
    double dt)
{
    State k1 = f(x, u);
    State k2 = f(x + 0.5 * dt * k1, u);
    State k3 = f(x + 0.5 * dt * k2, u);
    State k4 = f(x + dt * k3, u);

    return x + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}
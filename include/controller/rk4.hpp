#pragma once

#include <Eigen/Dense>

template <typename DynamicsFunction>
Eigen::VectorXd rk4Step(
    DynamicsFunction f,
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& t,
    double dt)
{
    State k1 = f(x, u);
    State k2 = f(x + 0.5 * dt * k1, t + 0.5 * dt);
    State k3 = f(x + 0.5 * dt * k2, t + 0.5 * dt);
    State k4 = f(x + dt * k3, t + dt);

    return x + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}
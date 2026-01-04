#pragma once
#include <Eigen/Dense>

struct linkInertia
{
    double mass;
    Eigen::Vector3d com;
    Eigen::Matrix3d inertia;
};
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "linearMpcHumanoid/general/generalizedFunctions.hpp"

//Computes the polynomial coefficientss for position, velocity and acceleration of the foot
//It assumes a 3th order polynomial for x and y trajectories
//A 7th order polynomial is assumed for the z trajectory

std::vector<Eigen::VectorXd> footCoeffTrajectory(const Eigen::Vector3d& currentPos,
    const Eigen::Vector3d& desPos,
    double stepHeight,
    double T);
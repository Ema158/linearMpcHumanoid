#pragma once
#include <vector>
#include <Eigen/Dense>

//Computes the polynomial coefficientss for position, velocity and acceleration of the foot
//It assumes a 3th order polynomial for x and y trajectories
//A 7th order polynomial is assumed for the z trajectory

std::vector<Eigen::VectorXd> footCoeffTrajectory(Eigen::VectorXd currentPos,
    Eigen::VectorXd desPos,
    double stepHeight,
    double T);
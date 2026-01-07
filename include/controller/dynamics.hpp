#pragma once
#include "controller/robotInfo.hpp"
#include "controller/generalizedFunctions.hpp"
#include "controller/linkInertia.hpp"
#include <Eigen/Dense>
#include <vector>

class dynamics{
public:

private:
    Eigen::MatrixXd spatialInertiaMatrix(linkInertia link); //compues the spatial inertia 6x6 matrix used in recursive algorithms
     
    std::vector<Eigen::MatrixXd> allSpatialInertiaMatrices(std::vector<linkInertia> links, robotInfo robot);
    Eigen::VectorXd computeC(robotInfo robot, std::vector<Eigen::MatrixXd> X, std::vector<Eigen::MatrixXd> I, bool isGravity); //..
    //..computes the Coriolles,centrigula,gravitational force vector.
    //it is the forward Newton-Euler with vD=0
};
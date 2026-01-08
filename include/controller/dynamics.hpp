#pragma once
#include "controller/robotInfo.hpp"
#include "controller/generalizedFunctions.hpp"
#include "controller/linkInertia.hpp"
#include <Eigen/Dense>
#include <vector>

class dynamics{
public:
    std::vector<Eigen::MatrixXd> allSpatialInertiaMatrices(robotInfo robot);
    Eigen::VectorXd computeC(robotInfo robot, std::vector<Eigen::MatrixXd> I, bool isGravity); //..
    Eigen::MatrixXd computeM(robotInfo robot, std::vector<Eigen::MatrixXd> I);//Composite rigid body algorithm
    Eigen::MatrixXd centroidalMatrix(Eigen::MatrixXd M, robotInfo robot); //computes the centroidal matrix of the centroidal model
private:
    Eigen::MatrixXd spatialInertiaMatrix(linkInertia link); //compues the spatial inertia 6x6 matrix used in recursive algorithms
     
    
    
    //..computes the Coriolles,centrigula,gravitational force vector.
    //it is the forward Newton-Euler with vD=0
};
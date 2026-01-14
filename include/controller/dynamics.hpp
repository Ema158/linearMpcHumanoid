#pragma once
#include "controller/Robot.hpp"
#include "controller/generalizedFunctions.hpp"
#include "controller/linkInertia.hpp"
#include <Eigen/Dense>
#include <vector>

class Dynamics{
    public:
    Dynamics();

    std::vector<Eigen::MatrixXd> allSpatialInertiaMatrices(
        const Robot& robot); //computes the spatial inertia 6x6 matrices of all frames

    void computeC(
        const Robot& robot,
        const std::vector<Eigen::MatrixXd>& I,
        bool isGravity); //..computes the Coriolles,centrigula,gravitational force vector.
        //it is the forward Newton-Euler with vD=0

    Eigen::MatrixXd computeM(
        const Robot& robot,
        const std::vector<Eigen::MatrixXd>& I);//Composite rigid body algorithm

    Eigen::MatrixXd centroidalMatrix(
        const Eigen::MatrixXd& M,
        const Robot& robot); //computes the centroidal matrix of the centroidal model
    
    Eigen::MatrixXd spatialInertiaMatrix(
        const linkInertia& link); //compues the spatial inertia 6x6 matrix of a fraame

    Eigen::VectorXd getC() const {return C_;}
    Eigen::MatrixXd getM() const {return M_;}
    Eigen::VectorXd getAG() const {return AG_;}
    
    private:
    Eigen::VectorXd C_;
    Eigen::MatrixXd M_;
    Eigen::MatrixXd AG_;
    Eigen::MatrixXd AGpqp_;
};
#pragma once
#include <Eigen/Dense>
#include "controller/Robot.hpp"
#include "controller/generalizedFunctions.hpp"
#include <cmath>

//Operational space is composed by 
// Position and orientation of each foot (12)
// Joint values for the arms and head (12)
// Orientation of the base in euler angles (3)
// Position of the center of mass (3)
namespace ik{
    Eigen::VectorXd desiredOperationalState(
        const Robot& robot,
        const Eigen::VectorXd& Rf,
        const Eigen::VectorXd& Lf,
        const Eigen::Vector3d& com); //

    void compute(
        Robot& robot,
        const Eigen::VectorXd& desOp);

    Eigen::VectorXd operationalState(
        const Robot& robot);

    Eigen::MatrixXd feetJacobian(
        const Robot& robot);

    Eigen::MatrixXd frameJacobian(
        const std::vector<Eigen::MatrixXd>& X,
        int frame,
        const Robot& robot);

    Eigen::MatrixXd comJacobian(
        const Robot& robot);

    Eigen::MatrixXd jacInvKinematics(
        const Robot& robot);
    
    Eigen::MatrixXd baseJacobian(
        const Eigen::Vector3d& pBase,
        const Eigen::Vector3d& pFrame);

    Eigen::Vector3d rotMatrixToEulerAngles(
        const Eigen::Matrix3d& R,
        const Eigen::Matrix3d& frame);//convertes a rotation matrix wrt world frame...
                                                                                    //to euler angles wrt world
                                                                                    //a reference frame is need to indicate the reference
                                                                                    //when eulerAnlges=0

    std::vector<Eigen::MatrixXd> allVelocityMatrices(
        const std::vector<Eigen::Matrix4d>& piTi);

    std::vector<Eigen::Matrix4d> parentTransMatrix(
        const std::vector<Eigen::Matrix4d>& T);
}
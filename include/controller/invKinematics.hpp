#pragma once
#include <Eigen/Dense>
#include "controller/robotInfo.hpp"
#include "controller/generalizedFunctions.hpp"
#include <cmath>

//Operational space is composed by 
// Position and orientation of each foot (12)
// Joint values for the arms and head (12)
// Orientation of the base in euler angles (3)
// Position of the center of mass (3)
class invKinematics{
public:
    Eigen::VectorXd desiredOperationalState(robotInfo robot, const Eigen::VectorXd Rf, const Eigen::VectorXd Lf, const Eigen::Vector3d com); //
    Eigen::VectorXd compute(robotInfo robot, Eigen::VectorXd desOp);
    Eigen::VectorXd operationalState(robotInfo robot);
    Eigen::MatrixXd feetJacobian(robotInfo robot);
    Eigen::MatrixXd frameJacobian(std::vector<Eigen::MatrixXd> X, int frame, robotInfo robot);
    Eigen::MatrixXd comJacobian(robotInfo robot);
    Eigen::MatrixXd jacInvKinematics(robotInfo robot);
private:
    
    Eigen::MatrixXd baseJacobian(Eigen::Vector3d pBase, Eigen::Vector3d pFrame);
    Eigen::Vector3d rotMatrixToEulerAngles(Eigen::Matrix3d R, Eigen::Matrix3d frame);//convertes a rotation matrix wrt world frame...
                                                                                    //to euler angles wrt world
                                                                                    //a reference frame is need to indicate the reference
                                                                                    //when eulerAnlges=0
};
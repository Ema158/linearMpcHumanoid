#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

Eigen::Matrix3d crossMatrix(Eigen::Vector3d v); //Convert a 3d vector into a skew symmetric matrix (cross product)

Eigen::MatrixXd velocityMatrix(Eigen::Matrix4d T); //given a transformation matrix that relates a frame wrt to its parent
                                                   //computes a 6x6 matrix that relates the spatial velocity of a body wrt to its parent

Eigen::Matrix4d inverseTransformationMatrix(Eigen::Matrix4d T);//Computes the inverse of a 4x4 transformation matrix  

Eigen::MatrixXd spatialCrossMatrix(Eigen::VectorXd v);//given a 6x1 spatial velocity
                                                     //computes the 6x6 matrix equal to the velocity spatial cross product

Eigen::MatrixXd spatialCrossMatrixForce(Eigen::VectorXd v);//given a 6x1 spatial velocity
                                                       //computes the 6x6 matrix equal to the force spatial cross product  
                                                       
Eigen::Matrix3d matrixAngularVelToEulerDot(Eigen::Vector3d eta); //Transformation matrix that transform angular velocity...
                                                              //to euler angles derivatives

Eigen::Matrix3d eulerAnglesToSO3(const Eigen::Vector3d& eulerAngles);     

Eigen::Vector3d rotMatToAxisAngle(const Eigen::Matrix3d& R);

Eigen::VectorXd findPolyCoeff(Eigen::MatrixXd Pos, Eigen::MatrixXd Vel, Eigen::MatrixXd Acc);

double polyval(const Eigen::VectorXd& poly, double x);

const Eigen::VectorXd polyder(const Eigen::VectorXd& poly);

void swapBaseVelocityAndRefToWorldFrame(const Eigen::MatrixXd& X01, Eigen::VectorXd& v);    
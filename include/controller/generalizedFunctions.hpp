#pragma once
#include <Eigen/Dense>

Eigen::Matrix3d crossMatrix(Eigen::Vector3d v); //Convert a 3d vector into a skew symmetric matrix (cross product)

Eigen::MatrixXd velocityMatrix(Eigen::Matrix4d T); //given a transformation matrix that relates a frame wrt to its parent
                                                   //computes a 6x6 matrix that relates the spatial velocity of a body wrt to its parent

Eigen::Matrix4d inverseTransformationMatrix(Eigen::Matrix4d T);//Computes the inverse of a 4x4 transformation matrix  

Eigen::MatrixXd spatialCrossMatrix(Eigen::VectorXd v);//given a 6x1 spatial velocity
                                                     //computes the 6x6 matrix equal to the spatial cross product
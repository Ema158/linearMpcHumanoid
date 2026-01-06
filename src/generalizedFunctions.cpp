#include "controller/generalizedFunctions.hpp"

Eigen::Matrix3d crossMatrix(Eigen::Vector3d v){
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero(3,3); 
    A <<    0, -v(2),  v(1),
         v(2),     0, -v(0), 
        -v(1),  v(0),     0;
    return A;
}

Eigen::MatrixXd velocityMatrix(Eigen::Matrix4d T){
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(6,6);
    Eigen::Matrix3d R = T.block(0,0,3,3).transpose(); //Rotation matrix
    Eigen::Vector3d p = T.block(0,3,3,1); //position vector
    X.block(0,0,3,3) = R;
    X.block(3,0,3,3) = -R*crossMatrix(p);
    X.block(3,3,3,3) = R;
    return X;
}

Eigen::Matrix4d inverseTransformationMatrix(Eigen::Matrix4d T){
    Eigen::Matrix4d Tinv = Eigen::Matrix4d::Zero(4,4);
    Tinv.block(0,0,3,3) = T.block(0,0,3,3).transpose();
    Tinv.block(0,3,3,1) = (-T.block(0,0,3,3).transpose())*T.block(0,3,3,1);
    Tinv(3,3) = 1;
    return Tinv;
}

Eigen::MatrixXd spatialCrossMatrix(Eigen::VectorXd v){
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(6,6);
    m.block(0,0,3,3) = crossMatrix(v.segment(0,3));
    m.block(3,0,3,3) = crossMatrix(v.segment(3,3));
    m.block(3,3,3,3) = crossMatrix(v.segment(0,3));
    return m;
}
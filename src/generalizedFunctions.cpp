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

Eigen::MatrixXd spatialCrossMatrixForce(Eigen::VectorXd v){
    Eigen::MatrixXd f = Eigen::MatrixXd::Zero(6,6);
    f = -spatialCrossMatrix(v).transpose();
    return f;
}

Eigen::Matrix3d matrixAngularVelToEulerDot(Eigen::Vector3d eta){
    //eta = [roll,pitch,yaw]
    Eigen::Matrix3d Omega = Eigen::Matrix3d::Zero(3,3);
    Omega << std::cos(eta(2))/std::cos(eta(1)),  std::sin(eta(2))/std::cos(eta(1)), 0,
                             -std::sin(eta(2)),                   std::cos(eta(2)), 0,
              std::cos(eta(2))*std::tan(eta(1)), std::sin(eta(2))*std::tan(eta(1)), 1;
    return Omega;
}

Eigen::Matrix3d eulerAnglesToSO3(const Eigen::Vector3d& rpy)
{
    const double roll  = rpy(0);
    const double pitch = rpy(1);
    const double yaw   = rpy(2);

    const double cr = std::cos(roll);
    const double sr = std::sin(roll);
    const double cp = std::cos(pitch);
    const double sp = std::sin(pitch);
    const double cy = std::cos(yaw);
    const double sy = std::sin(yaw);

    Eigen::Matrix3d R;

    R <<  cy * cp,  cy * sp * sr - sy * cr,  cy * sp * cr + sy * sr,
          sy * cp,  sy * sp * sr + cy * cr,  sy * sp * cr - cy * sr,
          -sp,      cp * sr,                 cp * cr;

    return R;
}
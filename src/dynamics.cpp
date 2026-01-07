#include <controller/dynamics.hpp>

Eigen::MatrixXd dynamics::spatialInertiaMatrix(linkInertia link){
    Eigen::MatrixXd I = Eigen::MatrixXd::Zero(6,6);
    I.block(0,0,3,3) = link.inertia - link.mass*crossMatrix(link.com)*crossMatrix(link.com);
    I.block(0,3,3,3) = link.mass*crossMatrix(link.com);
    I.block(3,0,3,3) = -link.mass*crossMatrix(link.com);
    I.block(3,3,3,3) = link.mass*Eigen::Matrix3d::Identity();
    return I;
}

std::vector<Eigen::MatrixXd> dynamics::allVelocityMatrices(robotInfo robot){
    std::vector<Eigen::Matrix4d> T = robot.getT();
    std::vector<Eigen::MatrixXd> X;
    std::vector<int> act = robot.actuatedFrames();
    X.resize(robot.getNumFrames());
    X[0] = velocityMatrix(T[0]); //Velocity matrix of the base frame
    for(int i=0;i<robot.getNumFrames();i++){
        if (act[i]!=0){ //Only compute the velocity transformation matrix of frames with actual joints
            X[i] = velocityMatrix(T[i]);
        }
    }
    return X;
}

std::vector<Eigen::MatrixXd> dynamics::allSpatialInertiaMatrices(std::vector<linkInertia> links, robotInfo robot){
    std::vector<Eigen::MatrixXd> I;
    std::vector<int> act = robot.actuatedFrames();
    I.resize(robot.getNumFrames());
    I[0] = spatialInertiaMatrix(robot.getLinks()[0]);
    for(int i=0;i<robot.getNumFrames();i++){
        if (act[i]!=0){ //Only compute the spatial inertia of frames with actual joints
            I[i] = spatialInertiaMatrix(links[i]);
        }
    }
    return I;
}

Eigen::VectorXd dynamics::computeC(robotInfo robot, std::vector<Eigen::MatrixXd> X, std::vector<Eigen::MatrixXd> I){
   Eigen::VectorXd C = Eigen::VectorXd::Zero(robot.getNumJoints());
   return C; 
}
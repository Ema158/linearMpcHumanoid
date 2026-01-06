#include "controller/invKinematics.hpp"

Eigen::VectorXd invKinematics::desiredOperationalState(robotInfo robot, const Eigen::VectorXd Rf, const Eigen::VectorXd Lf, const Eigen::Vector3d com){
    Eigen::VectorXd Qd = Eigen::VectorXd::Zero(robot.getNumJoints());
    Eigen::VectorXd q = robot.getJoints();
    Qd.segment(0,6) = Rf;
    Qd.segment(6,6) = Lf;
    Qd.segment(6,12) = q.segment(12,12);
    return Qd;
}

Eigen::VectorXd compute(robotInfo robot, Eigen::VectorXd desOp){
    Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.getNumJoints());
    return q;
}

Eigen::VectorXd operationalState(robotInfo robot){
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(robot.getNumJoints());
    std::vector<Eigen::Matrix4d> T = robot.getT();
    Eigen::VectorXd q = robot.getJoints(); 
    Q.segment(0,3) = T[7].block(0,3,3,1); //Position of the right foot
    //The frame attached to the foot does not coincide when the inertial frame when q=0 (rotRF(q=0) =! identity)
    //Then, if we want Roll,Pitch,Yaw angles of the foot wrt inertial frame, the rotation matrix of the foot...
    //... must change such that when q=0 rotRF(q=0) = identity
    Eigen::Matrix3d rotRF; //Rotation matrix of the right foot such that when q=0 rotRF(q=0) = identity
    rotRF = T[7].block(0,0,3,3)*robot.Rf_q0;
    Q(5) = std::atan2(rotRF(1,0),rotRF(0,0)); //Yaw Right foot
    Q(4) = std::atan2(-rotRF(2,0),std::cos(Q(5))*rotRF(0,0)+std::sin(Q(5))*rotRF(1,0)); //Pitch right foot
    Q(3) = std::atan2(std::sin(Q(5))*rotRF(0,2)-std::cos(Q(5))*rotRF(1,2),-std::sin(Q(5))*rotRF(0,1)+std::cos(Q(5))*rotRF(1,1)); //Roll right foot

    Q.segment(6,3) = T[14].block(0,3,3,1); //Position of the left foot
    Eigen::Matrix3d rotLF; //Rotation matrix of the left foot such that when q=0 rotRF(q=0) = identity
    rotLF = T[14].block(0,0,3,3)*robot.Lf_q0;
    Q(11) = std::atan2(rotLF(1,0),rotLF(0,0)); //Yaw Right foot
    Q(10) = std::atan2(-rotLF(2,0),std::cos(Q(11))*rotLF(0,0)+std::sin(Q(11))*rotLF(1,0)); //Pitch right foot
    Q(9) = std::atan2(std::sin(Q(11))*rotLF(0,2)-std::cos(Q(11))*rotLF(1,2),-std::sin(Q(11))*rotLF(0,1)+std::cos(Q(11))*rotLF(1,1));
    
    Q.segment(12,12) = q.segment(12,robot.getNumJoints()); //arms and head joints
    Q.segment(24,3) = q.segment(3,3);
    //Q.segment(27,3) = robot.getCoM();
    return Q;
}
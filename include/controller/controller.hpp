#pragma once
#include <Eigen/Dense>
#include "controller/Robot.hpp"
#include "controller/mpcLinearPendulum.hpp"
#include "controller/zmpGeneration.hpp"
#include "controller/Task.hpp"
#include "controller/Dynamics.hpp"
#include "controller/Clock.hpp"

/*
state = [q, v]
q = [0p,0eta,qJ]
    0p ->   base position wrt world frame
    0eta -> base orientation wrt world frame (Euler angles)
    qJ ->   actual joints

v = [0v,0w, qDJ]
    0v ->  linear velocity of the base
    0w ->  angular velocity of the base
    qDJ -> actual joints velocity
*/

class Controller {
public:
    Controller(Robot& robot);

    void stand();

    void computeComVelocity(Eigen::VectorXd v); //Uses state and centroidal matrix to compute the velocity of the center of mass
                               //Other option is to use the center of mass Jacobian of invKinematics

private:
    double simulationTime_ = 4;
    double t_ = 0;
    double dt_ = 0.01;
    Eigen::VectorXd tau_; //24 dimention vector of torques
    Eigen::VectorXd state_; // 60 dimention vector of the current configuration and velocity
    Robot& robot_;
    Dynamics dyn_;
    Eigen::Vector3d comVel_ = Eigen::Vector3d::Zero();
};
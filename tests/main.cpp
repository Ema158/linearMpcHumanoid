#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/Robot.hpp"
#include "controller/invKinematics.hpp"
#include "controller/dynamics.hpp"
#include "controller/mpcLinearPendulum.hpp"
#include <iostream>

int main() {
    Robot nao;
    Eigen::VectorXd Rf = Eigen::VectorXd::Zero(6);
    Rf(1) = -0.05;
    Eigen::VectorXd Lf = Eigen::VectorXd::Zero(6);
    Lf(1) = 0.05;
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    com(2) = 0.26;
    Eigen::VectorXd desOp = ik::desiredOperationalState(nao,Rf,Lf,com);
    ik::compute(nao, desOp);
    Mpc3dLip mpc;

    
    return 0;
}
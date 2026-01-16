#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/Robot.hpp"
#include "controller/invKinematics.hpp"
#include "controller/Dynamics.hpp"
#include "controller/mpcLinearPendulum.hpp"
#include "controller/zmpGeneration.hpp"
#include "controller/Clock.hpp"
#include <iostream>

int main() {
    Robot nao;
    Kinematics ik;
    Clock clock;
    Eigen::VectorXd Rf = Eigen::VectorXd::Zero(6);
    Rf(1) = -0.05;
    Eigen::VectorXd Lf = Eigen::VectorXd::Zero(6);
    Lf(1) = 0.05;
    Eigen::Vector3d com = Eigen::Vector3d::Zero();

    com(2) = 0.26;
    com(1) = 0.01;
    com(0) = -0.005;
    Eigen::VectorXd desOp = ik.desiredOperationalState(nao,Rf,Lf,com);
    ik.compute(nao, desOp);
    
    double timeHorizon = 0.5;
    Mpc3dLip mpc(clock.getTimeStep(), timeHorizon, nao.getCoM()(2));
    Controller controller(nao,mpc);
    controller.stand();

    return 0;
}
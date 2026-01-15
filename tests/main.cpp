#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/Robot.hpp"
#include "controller/invKinematics.hpp"
#include "controller/Dynamics.hpp"
#include "controller/mpcLinearPendulum.hpp"
#include "controller/zmpGeneration.hpp"
#include <iostream>

int main() {
    Robot nao;
    Kinematics ik;
    
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
    
    Controller controller(nao);
    controller.stand();
    //test mpc
    /*Eigen::VectorXd newState = Eigen::VectorXd(6);
    Eigen::Vector2d vel_xy = Eigen::Vector2d::Zero();
    for (int i=0;i<100;i++){
        newState = mpc.compute(com_xy,vel_xy,zmp.getZmpXRef(),zmp.getZmpYRef());
        com_xy << newState(0),newState(3);
        vel_xy << newState(1),newState(4);
        std::cout<<newState<<std::endl<<std::endl;
    }*/
    return 0;
}
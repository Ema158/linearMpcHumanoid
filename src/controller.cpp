#include "controller/controller.hpp"
#include <iostream>
#include <cmath>

Controller::Controller(
    Robot& robot)
    :
    robot_(robot) 
{
    std::cout<<"Controller Initiated" << std::endl;
    std::cout<<"Initial conditions of the center of mass: " << std::endl << "c = " << robot_.getCoM();

}

void Controller::stand() 
{
     Mpc3dLip mpc;
     ZMP zmp(Task::Stand,simulationTime_,dt_,SupportFoot::Double);
     Eigen::Vector2d com_xy;
     com_xy << robot_.getCoM()(0),robot_.getCoM()(1);
     Eigen::Vector2d vel_xy = Eigen::Vector2d::Zero(2);
     Eigen::VectorXd newState = Eigen::VectorXd::Zero(6);
     while(std::abs(t_ - simulationTime_) > 0.01)
     {
        newState = mpc.compute(com_xy,vel_xy,zmp.getZmpXRef(),zmp.getZmpYRef());
        com_xy << newState(0),newState(3);
        vel_xy << newState(1),newState(4);
        std::cout<<newState<<std::endl<<std::endl;

        t_ = t_ + dt_;
     }   
}



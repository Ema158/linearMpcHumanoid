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
    //Initial State
    state_.resize(2*robot_.getNumJoints());
    state_.segment(0,robot_.getNumJoints()) = robot_.getJoints();
    state_.segment(robot_.getNumJoints(),robot_.getNumJoints()) = robot_.getJointsVelocity();
    //Initial torques
    tau_ = Eigen::VectorXd::Zero(robot_.getNumJoints());
}

void Controller::stand() 
{
     Clock clock;
     double timeHorizon = 0.5;
     Mpc3dLip mpc(clock.getTimeStep(), timeHorizon, robot_.getCoM()(2));
     ZMP zmp(Task::Stand,simulationTime_,dt_,SupportFoot::Double);
     Eigen::Vector2d com_xy;
     Eigen::Vector2d comVel_xy;
     Eigen::Vector2d vel_xy = Eigen::Vector2d::Zero(2);
     Eigen::VectorXd newState = Eigen::VectorXd::Zero(6);
     com_xy << robot_.getCoM()(0),robot_.getCoM()(1);
     while(std::abs(clock.getTime() - simulationTime_) > 0.01)
     {
        dyn_.computeAll(robot_); //Computes M,C,AG,AGpqp

        //com_xy << robot_.getCoM()(0),robot_.getCoM()(1); //Current position center of mass
        computeComVelocity(state_.segment(robot_.getNumJoints(),robot_.getNumJoints())); //Current velocity center of mass
        //comVel_xy << comVel_(0),comVel_(1);
        
        newState = mpc.compute(com_xy,comVel_xy,zmp.getZmpXRef(),zmp.getZmpYRef(), clock.getTime());
        com_xy << newState(0),newState(3);
        comVel_xy << newState(1),newState(4);
        std::cout<<newState<<std::endl<<std::endl;

        clock.step();
     }   
}

void Controller::computeComVelocity(Eigen::VectorXd v)
{
    // The order of the velocity variables in the centroidal matrix is [angular vel, linear velocity]
    // v is in the order [linear velocity, angular velocity]. We need to change
    Eigen::Vector3d temp = v.segment(0,3);
    v.segment(0,3) = v.segment(3,3);
    v.segment(3,3) = temp;

    //Also base velocity in v is wrt to world frame
    //In centroidal matrix base velocity in v is wrt base frame
    v.segment(0,6) = robot_.getX()[0]*v.segment(0,6);

    comVel_ = dyn_.getAG().block(3,0,3,robot_.getNumJoints())*v; //AG*v returns linear and angular momentum
                                                                //We just need linear

    comVel_ = comVel_/robot_.getMass(); //We divided by mass to obtain velocity                                                                
}



#include "simulators/MujocoSim.hpp"
#include "simulators/MujocoViewer.hpp"

#include "linearMpcHumanoid/robotInfo/Robot.hpp"
#include "linearMpcHumanoid/controller/controller.hpp"
#include "linearMpcHumanoid/controller/invKinematics.hpp"
#include "linearMpcHumanoid/controller/mpcLinearPendulum.hpp"
#include "linearMpcHumanoid/general/Clock.hpp"
#include "linearMpcHumanoid/general/Task.hpp"

#include <iostream>
#include <Eigen/Dense>

void stand(Robot& robot, Controller& controller, Clock& clock, MujocoSim& sim);

Eigen::VectorXd dynamics(const Eigen::VectorXd& state, double t, Robot& Robot, Controller& controller);

Eigen::MatrixXd relabelMujocoMatrix(Robot& robot);

int main() {
  double simulationTime = 0.05;
  double timeStep = 0.01;

  //Desired initial configuration for the simulation
  Robot nao;
  Kinematics ik;
  Clock clock(timeStep,simulationTime);

  //ZMP trajectory for a stand task (in the center of the support zone for all time)
  ZMP zmp(Task::Stand,simulationTime,timeStep,SupportFoot::Double);
    
  //Initial position of the feet for simulation
  Eigen::VectorXd Rf = Eigen::VectorXd::Zero(6);
  Rf(1) = -0.05;
  Eigen::VectorXd Lf = Eigen::VectorXd::Zero(6);
  Lf(1) = 0.05;
    
  //Initial position of the center of mass for simulation
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  com << 0.02, 0.00, 0.26 ;
    
  //Inverse kinematics to compute the initial joint configuration
  Eigen::VectorXd desOp = ik.desiredOperationalState(nao,Rf,Lf,com);
  ik.compute(nao, desOp);
  Eigen::VectorXd q0 = nao.getJoints();

  //Prediction time of the model predictive controller of the linear inverted pendulum model and initialization
  double timeHorizon = 0.5;
  Mpc3dLip mpc(clock.getTimeStep(), timeHorizon, nao.getCoM()(2));

  //Desired trajectory for the feet during simulation
  //No movement for stand
  std::vector<Eigen::VectorXd> rFCoeff;
  std::vector<Eigen::VectorXd> lFCoeff;
  Eigen::Vector3d currentPos;
  Eigen::Vector3d desPos;

  currentPos << 0, -0.05, 0;
  desPos << 0, -0.05, 0;
  double stepHeight = 0;
  rFCoeff = footCoeffTrajectory(currentPos, desPos, stepHeight, simulationTime);

  currentPos << 0, 0.05, 0;
  desPos << 0, 0.05, 0;
  stepHeight = 0;
  lFCoeff = footCoeffTrajectory(currentPos, desPos, stepHeight, simulationTime);
    
  Controller controller(nao,mpc,zmp,rFCoeff,lFCoeff);
  
  //Create simulation and viewer
  MujocoSim sim("../models/nao.xml", q0);
  MujocoViewer viewer(sim);

  Eigen::VectorXd q_test(24);
  Eigen::MatrixXd L = relabelMujocoMatrix(nao);

  q_test = L*q0.segment(6,nao.getNumActualJoints());           

  //Simulation
  viewer.run([&]() {
        stand(nao,controller,clock,sim);
        //q0 = nao.getJoints();
        q_test = L*q0.segment(6,nao.getNumActualJoints()); 
        sim.applyJointPositions(q_test);
    }, clock);

  /*Eigen::VectorXd tau_test(24);
  Eigen::VectorXd tau0(24);
  viewer.run([&]() {
        stand(nao,controller,clock);
        tau0 = controller.getTorques();
        tau_test << tau0(22), tau0(23), 
            tau0(6), tau0(7), tau0(8), tau0(9), tau0(10), tau0(11),
            tau0(0), tau0(1), tau0(2), tau0(3), tau0(4), tau0(5),
            -tau0(17), tau0(18), tau0(19), tau0(20), tau0(21),
            tau0(12), tau0(13), tau0(14), tau0(15), tau0(16); 
        sim.applyJointPositions(tau_test);
    }, clock);*/
 
  return 0;
}

void stand(Robot& robot, Controller& controller, Clock& clock, MujocoSim& sim) 
{
    int n = robot.getNumJoints();
    Eigen::VectorXd state(2*n);

    Eigen::VectorXd currentPos(n+1);
    Eigen::VectorXd currentVel(n);
    
    //Get state from rk4 integration (no feedback)
    state.segment(0,n) = robot.getJoints();
    state.segment(n,n) = robot.getJointsVelocity();

    //Get state from mujoco data (feedback), no orientation for now
    //sim.getMujocoState(currentPos, currentVel);
    //state.segment(0,3) = currentPos.segment(0,3);
    //state(2) -= 0.035;
    //Eigen::VectorXd relabelJoints(24); 
    //relabelJoints = (relabelMujocoMatrix(robot).transpose()) * (currentPos.segment(7,24));
    //state.segment(6,6) = relabelJoints.segment(0,6);
    //state.segment(n,3) = currentVel.segment(0,3);
  
    state = rk4Step(
        [&](const Eigen::VectorXd& x, double t)
        {
            return dynamics(x, t, robot, controller);
        },
    state,
        clock.getTime(),
        clock.getTimeStep()
        );

    //std::cout<<robot.getCoM()(1)<<std::endl<<std::endl;
    //std::cout << "Mujoco x position = " << currentPos(0) << " RK4 position = " << state(0) << std::endl;
    //std::cout << "Mujoco y position = " << currentPos(1) << " RK4 position = " << state(1) << std::endl;
    //std::cout << "Mujoco z position = " << currentPos(2)- 0.035 << " RK4 position = " << state(2) << std::endl;
    //for(int i=0;i<24;i++){
     //   std::cout << "Mujoco q"<< i+1 << " = " << currentPos(i+6) << " RK4 q"<<i+1 <<"= " << state(i+6) << std::endl;
        //std::cout << "Mujoco q"<< i+1 << " = " << relabelJoints(i) << " RK4 q"<<i+1 <<"= " << state(i+6) << std::endl;
    //}
    std::cout<<std::endl;
    /*std::cout << "Mujoco x velocity = " << currentVel(0) << " RK4 position = " << state(n) << std::endl;
    std::cout << "Mujoco y position = " << currentVel(1) << " RK4 position = " << state(n+1) << std::endl;
    std::cout << "Mujoco z position = " << currentVel(2) << " RK4 position = " << state(n+2) << std::endl << std::endl;*/

    clock.step();    
}

Eigen::VectorXd dynamics(const Eigen::VectorXd& state, double t, Robot& robot, Controller& controller)
{
    const int n = robot.getNumJoints();

    Eigen::VectorXd q = state.segment(0,n);
    Eigen::VectorXd qD = state.segment(n,n);

    ControllerInput in;
    in.q    = state.segment(0,n);
    in.dq   = state.segment(n,n);
    in.time = t;

    controller.standStep(in);

    WBCOutput out = controller.WBC(state, t);

    //The linear and angular velocities in the state vector are spatial velocities
    //Before integrating linear velocity must change to classical definition
    //Definition of 0v1: is the linear velocity of a point joined to body 1 that currently coincide with the origin of frame 0
    //We transform this linear velocity with the classical lineal velocity transformation
    qD.segment(0,3) += crossMatrix(qD.segment(3,3)) * q.segment(0,3); // 0v1 (classic) = 0v1(spatial) + 0w1 X 0p1
    // Classical and spatial angular velocity are the same

    //The rotation of the base is represented as Euler Angles
    //We need to transform angular velocity to Euler Angles rate before integrating
    qD.segment(3,3) =  matrixAngularVelToEulerDot(q.segment(3,3))*qD.segment(3,3); // 0eta1 = Omega*0w1

    Eigen::VectorXd xp(2 * n);
    xp.head(n) = qD;// q̇
    xp.tail(n) = out.qpp;// q̈
    return xp;
}

Eigen::MatrixXd relabelMujocoMatrix(Robot& robot)
{
    Eigen::MatrixXd L(robot.getNumActualJoints(), robot.getNumActualJoints());
    int numHeadJoints = 2;
    int numLegJoints = 6;
    int numArmJoints = 5;

    L.block(0,22,numHeadJoints,numHeadJoints) = Eigen::Matrix2d::Identity();

    L.block(2,6,numLegJoints,numLegJoints) = Eigen::MatrixXd::Identity(numLegJoints,numLegJoints);

    L.block(8,0,numLegJoints,numLegJoints) = Eigen::MatrixXd::Identity(numLegJoints,numLegJoints);

    L.block(14,17, numArmJoints, numArmJoints) = Eigen::MatrixXd::Identity(numArmJoints,numArmJoints);
    L(14,17) = -1;

    L.block(19,12, numArmJoints, numArmJoints) = Eigen::MatrixXd::Identity(numArmJoints,numArmJoints);
    return L;
}
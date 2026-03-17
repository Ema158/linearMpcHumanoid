#include "simulators/MujocoSim.hpp"
#include "simulators/MujocoViewer.hpp"
#include "simulators/MujocoDataToControllerData.hpp"

#include "linearMpcHumanoid/robotInfo/Robot.hpp"
#include "linearMpcHumanoid/controller/controller.hpp"
#include "linearMpcHumanoid/controller/invKinematics.hpp"
#include "linearMpcHumanoid/controller/mpcLinearPendulum.hpp"
#include "linearMpcHumanoid/general/Clock.hpp"
#include "linearMpcHumanoid/general/Task.hpp"
#include "linearMpcHumanoid/general/generalizedFunctions.hpp"
#include "linearMpcHumanoid/controller/ContactState.hpp"

#include <iostream>
#include <Eigen/Dense>

#include <chrono>

int main() {
  double simulationTime = 0.5;
  double timeStep = 0.01;
  ContactState contact(Contact::Both);
  Task task = Task::Walk;

  //Desired initial configuration for the simulation
  Robot nao;
  Kinematics ik;
  Clock clock(timeStep,simulationTime);
    
  //Initial position of the feet for simulation
  Eigen::VectorXd Rf = Eigen::VectorXd::Zero(6);
  Rf(1) = -0.05;
  Rf(2) = 0.00;
  Eigen::VectorXd Lf = Eigen::VectorXd::Zero(6);
  Lf(1) = 0.05;
  Lf(2) = 0.00;
    
  //Initial position of the center of mass for simulation
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  com << 0.00, 0.00, 0.26 ;
    
  //Inverse kinematics to compute the initial joint configuration
  Eigen::VectorXd desOp = ik.desiredOperationalState(nao,Rf,Lf,com);
  ik.compute(nao, desOp);
  Eigen::VectorXd q0 = nao.getJoints();

  //Prediction time of the model predictive controller of the linear inverted pendulum model and initialization
  double timeHorizon = 0.5;
  Mpc3dLip mpc(clock.getTimeStep(), timeHorizon, nao.getCoM()(2));

  //ZMP trajectory for a stand task (in the center of the support zone for all time)
  //ZMP zmp(task,simulationTime,timeStep,contact.get()); //Balance
  
  GaitParameters gaitParameters;
  ZMP zmp(task, timeStep, gaitParameters, contact.get()); //Walk
  
  //Desired trajectory for the feet during simulation
  //No movement for stand
  std::vector<Eigen::VectorXd> rFCoeff;
  std::vector<Eigen::VectorXd> lFCoeff;
  Eigen::Vector3d currentPos;
  Eigen::Vector3d desPos;

  currentPos << Rf(0), Rf(1), Rf(2);
  desPos << 0, -0.05, 0.00;
  double stepHeight = 0.00;
  rFCoeff = footCoeffTrajectory(currentPos, desPos, stepHeight, simulationTime);

  currentPos << Lf(0), Lf(1), Lf(2);
  desPos << 0, 0.05, 0.0;
  stepHeight = 0.0;
  lFCoeff = footCoeffTrajectory(currentPos, desPos, stepHeight, simulationTime);
    
  Controller controller(nao,mpc,zmp,rFCoeff,lFCoeff,contact);
  
  //Create simulation and viewer
  MujocoSim sim("../models/nao.xml", q0);
  MujocoViewer viewer(sim);       

  viewer.run([&]() {
        auto start = std::chrono::high_resolution_clock::now();

        //updateController(nao, controller, clock, sim);
        updateWalkController(nao, controller, clock, sim, zmp);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Elapsed time: " << elapsed.count() << " seconds\n"; 
    }, clock);
 
  return 0;
}






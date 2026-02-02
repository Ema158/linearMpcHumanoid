#include "simulators/MujocoSim.hpp"
#include "simulators/MujocoViewer.hpp"

#include "linearMpcHumanoid/robotInfo/Robot.hpp"
#include "linearMpcHumanoid/controller/controller.hpp"
#include "linearMpcHumanoid/controller/invKinematics.hpp"
//#include "linearMpcHumanoid/controlller/Dynamics.hpp"
#include "linearMpcHumanoid/controller/mpcLinearPendulum.hpp"
#include "linearMpcHumanoid/general/Clock.hpp"
#include "linearMpcHumanoid/general/Task.hpp"

#include <iostream>
#include <Eigen/Dense>

void stand(Robot& robot, Controller& controller, Clock& clock);

Eigen::VectorXd dynamics(const Eigen::VectorXd& state, double t, Robot& Robot, Controller& controller);

int main() {
  double simulationTime = 4;
  double timeStep = 0.01;
  Clock clock(timeStep,simulationTime);
  //Desired initial configuration for the simulation
  Robot nao;
  Kinematics ik;

  //ZMP trajectory for a stand task (in the center of the support zone for all time)
  ZMP zmp(Task::Stand,simulationTime,timeStep,SupportFoot::Double);
    
  //Initial position of the feet for simulation
  Eigen::VectorXd Rf = Eigen::VectorXd::Zero(6);
  Rf(1) = -0.05;
  Eigen::VectorXd Lf = Eigen::VectorXd::Zero(6);
  Lf(1) = 0.05;
    
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
  q_test << q0(28), q0(29), 
            q0(12), q0(13), q0(14), q0(15), q0(16), q0(17),
            q0(6), q0(7), q0(8), q0(9), q0(10), q0(11),
            -q0(23), q0(24), q0(25), q0(26), q0(27),
            q0(18), q0(19), q0(20), q0(21), q0(22);            

  //Simulation
  viewer.run([&]() {
        stand(nao,controller,clock);
        q0 = nao.getJoints();
        q_test << q0(28), q0(29), 
            q0(12), q0(13), q0(14), q0(15), q0(16), q0(17),
            q0(6), q0(7), q0(8), q0(9), q0(10), q0(11),
            -q0(23), q0(24), q0(25), q0(26), q0(27),
            q0(18), q0(19), q0(20), q0(21), q0(22); 
        sim.applyJointPositions(q_test);
    }, clock);

  return 0;
}

void stand(Robot& robot, Controller& controller, Clock& clock) 
{
    int n = robot.getNumJoints();
    Eigen::VectorXd state(2*n);
    state.segment(0,n) = robot.getJoints();
    state.segment(n,n) = robot.getJointsVelocity();
    Dynamics dyn;
  
    state = rk4Step(
        [&](const Eigen::VectorXd& x, double t)
        {
            return dynamics(x, t, robot, controller);
        },
    state,
        clock.getTime(),
        clock.getTimeStep()
        );

    std::cout<<state(0)<<std::endl<<std::endl;
    robot.updateState(state.segment(0,n));
    //robot.updateStateVelocity(state.segment(n,n), dyn.getAG());
    clock.step();    
}

Eigen::VectorXd dynamics(const Eigen::VectorXd& state, double t, Robot& robot, Controller& controller)
{
    const int n = robot.getNumJoints();

    Eigen::VectorXd q = state.segment(0,n);
    Eigen::VectorXd qD = state.segment(n,n);

    ControllerInput in;
    in.q    = robot.getJoints();
    in.dq   = robot.getJointsVelocity();
    in.time = t;

    Eigen::VectorXd tau = controller.standStep(in);

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
#include "simulators/MujocoSim.hpp"
#include "simulators/MujocoViewer.hpp"

#include "linearMpcHumanoid/robotInfo/Robot.hpp"
#include "linearMpcHumanoid/controller/controller.hpp"
#include "linearMpcHumanoid/controller/invKinematics.hpp"
#include "linearMpcHumanoid/controller/mpcLinearPendulum.hpp"
#include "linearMpcHumanoid/general/Clock.hpp"
#include "linearMpcHumanoid/general/Task.hpp"
#include "linearMpcHumanoid/general/generalizedFunctions.hpp"

#include <iostream>
#include <Eigen/Dense>

void stand(Robot& robot, Controller& controller, Clock& clock, MujocoSim& sim);

Eigen::VectorXd dynamics(const Eigen::VectorXd& state, double t, Robot& Robot, Controller& controller);

Eigen::MatrixXd relabelMujocoMatrix(Robot& robot);

Eigen::MatrixXd fromMujocoFrameToControllerFrame(const Eigen::MatrixXd& R01,const Eigen::Vector3d& p01);

int main() {
  double simulationTime = 5;
  double timeStep = 0.005;

  //Desired initial configuration for the simulation
  Robot nao;
  Kinematics ik;
  Clock clock(timeStep,simulationTime);

  //ZMP trajectory for a stand task (in the center of the support zone for all time)
  ZMP zmp(Task::Stand,simulationTime,timeStep,SupportFoot::Double);
    
  //Initial position of the feet for simulation
  Eigen::VectorXd Rf = Eigen::VectorXd::Zero(6);
  Rf(1) = -0.05;
  Rf(2) = 0.00;
  Eigen::VectorXd Lf = Eigen::VectorXd::Zero(6);
  Lf(1) = 0.05;
  Lf(2) = 0.00;
    
  //Initial position of the center of mass for simulation
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  com << 0.00, 0.02, 0.26 ;
    
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

  Eigen::VectorXd v0(24);
  Eigen::VectorXd v_test(24);

  Eigen::VectorXd tau_test = Eigen::VectorXd::Zero(24);
  Eigen::VectorXd tau0(24);
  viewer.run([&]() {
        stand(nao,controller,clock,sim);

        q0 = nao.getJoints();
        tau0 = controller.getTorques();

        tau_test = L*tau0; 

        std::vector<jointsIndex> joints = {
            jointsIndex::HeadYaw,
            jointsIndex::HeadPitch,
            jointsIndex::LHipYawPitch,
            jointsIndex::LHipRoll,
            jointsIndex::LHipPitch,
            jointsIndex::LKneePitch,
            jointsIndex::LAnklePitch,
            jointsIndex::LAnkleRoll,
            jointsIndex::RHipYawPitch,
            jointsIndex::RHipRoll,
            jointsIndex::RHipPitch,
            jointsIndex::RKneePitch,
            jointsIndex::RAnklePitch,
            jointsIndex::RAnkleRoll,
            jointsIndex::LShoulderPitch,
            jointsIndex::LShoulderRoll,
            jointsIndex::LElbowYaw,
            jointsIndex::LElbowRoll,
            jointsIndex::LWristYaw,
            jointsIndex::RShoulderPitch,
            jointsIndex::RShoulderRoll,
            jointsIndex::RElbowYaw,
            jointsIndex::RElbowRoll,
            jointsIndex::RWristYaw,
        };
        sim.applyTorquesV2(joints, tau_test);
    }, clock);
 
  return 0;
}

void stand(Robot& robot, Controller& controller, Clock& clock, MujocoSim& sim) 
{
    Dynamics dyn;
    Kinematics kin;
    int n = robot.getNumJoints();
    Eigen::VectorXd state(2*n);

    Eigen::VectorXd currentPos(n+1);
    Eigen::VectorXd currentVel(n);
    
    //Get state from rk4 integration (no feedback)
    state.segment(0,n) = robot.getJoints();
    state.segment(n,n) = robot.getJointsVelocity();

    //Get state from mujoco data (feedback)
    sim.getMujocoState(currentPos, currentVel);
    
    //Orientation in mujoco is represented as a quaternion
    //My current controller need euler angles
    Eigen::VectorXd quat = currentPos.segment(3,4);
    Eigen::Vector3d eulerAngles = quaternionToEulerAngles(quat);

    //===================================Base feedback
    //Base frame of mujoco and my controller are not the same
    //We need to change base position of mujoco to base position of my controller
    
    Eigen::Matrix3d R01 = quaternionToRotationMatrix(quat);
 
    Eigen::MatrixXd T01 = fromMujocoFrameToControllerFrame(R01, currentPos.segment(0,3));

    state.segment(0,3) = T01.block(0,3,3,1);
    state.segment(3,3) = eulerAngles;

    //==================================Base velocity feedback
    //Mujoco returns the velocity (linear and angular) of the base wrt world frame using classic definition of velocity
    //My controller need the spatial definition of velocity
    
    // 1v1 (spatial) = 1v1 (classic), same for angular
    // The first 1v1 (classic) = 1R0 0v1, 1w1 (classic) = 1R0 0w1 (classic)
    Eigen::VectorXd baseVelocityBaseFrame = Eigen::VectorXd::Zero(6);
    baseVelocityBaseFrame.segment(0,3) = R01.transpose()*currentVel.segment(3,3); //angular velocity first
    baseVelocityBaseFrame.segment(3,3) = R01.transpose()*currentVel.segment(0,3); //linear velocity

    //Now we use the transformation matrix for spatial velocity 
    // [0w1,0v1] (spatial) = 0X1 [1w1,1v1] (spatial)
    Eigen::VectorXd baseVelocityWorldFrameSpatial = Eigen::VectorXd::Zero(6);
    baseVelocityWorldFrameSpatial = (robot.getX()[0].inverse()) * baseVelocityBaseFrame;

    state.segment(n,3) = baseVelocityWorldFrameSpatial.segment(3,3); //linear
    state.segment(n+3,3) = baseVelocityWorldFrameSpatial.segment(0,3); //angular

    //=====================================Joints feedback
    Eigen::VectorXd relabelJoints(24); 
    assert(currentPos.size() >= 31);
    
    relabelJoints = (relabelMujocoMatrix(robot).transpose()) * (currentPos.segment(7,24));
    state.segment(6,24) = relabelJoints;

    Eigen::VectorXd relabelJointsVel(24); 
    relabelJointsVel = (relabelMujocoMatrix(robot).transpose()) * (currentVel.segment(6,24));
    state.segment(n+6,24) = relabelJointsVel;

    std::cout << "t = " << clock.getTime() << std::endl;
    std::cout << "Mujoco x com = " << sim.getCoM()(0) << " Controller x com = " << robot.getCoM()(0) << std::endl;
    std::cout << "Mujoco y com = " << sim.getCoM()(1) << " Controller y com = " << robot.getCoM()(1) << std::endl;
    std::cout << "Mujoco z com = " << sim.getCoM()(2) << " Controller z com = " << robot.getCoM()(2) << std::endl;  
    std::cout<< std::endl;

    ControllerInput in;
    in.q    = state.segment(0,n);
    in.dq   = state.segment(n,n);
    in.time = clock.getTime();
    controller.standStep(in);

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

    WBCOutput out = controller.WBC(t);

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
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(robot.getNumActualJoints(), robot.getNumActualJoints());
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

Eigen::MatrixXd fromMujocoFrameToControllerFrame(const Eigen::MatrixXd& R01,const Eigen::Vector3d& p01)
{
    Eigen::MatrixXd T01Mujoco = Eigen::MatrixXd::Zero(4,4); //Transformation matrix mujoco base frame wrt world frame
    T01Mujoco.block(0,0,3,3) = R01;
    T01Mujoco.block(0,3,3,1) = p01;
    T01Mujoco(3,3) = 1;

    Eigen::MatrixXd TMujoco_Controller = Eigen::MatrixXd::Zero(4,4); //Transformation matrix controller base frame wrt mujoco base frame
    TMujoco_Controller.block(0,0,3,3) = Eigen::Matrix3d::Identity();
    TMujoco_Controller(0,3) = 0;
    TMujoco_Controller(1,3) = 0;
    TMujoco_Controller(2,3) = -0.035;
    TMujoco_Controller(3,3) = 1;

    Eigen::MatrixXd T01Controller = T01Mujoco*TMujoco_Controller; //Transformation matrix controller base frame using mujoco data
    return T01Controller;
}
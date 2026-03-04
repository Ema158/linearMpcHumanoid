#include "simulators/MujocoDataToControllerData.hpp"

Eigen::VectorXd updateState(const Robot& robot, Clock& clock, MujocoSim& sim) 
{
    int n = robot.getNumJoints();
    Eigen::VectorXd state(2*n);

    Eigen::VectorXd currentPos(n+1);
    Eigen::VectorXd currentVel(n);

    //Get state from mujoco data (feedback)
    sim.getMujocoState(currentPos, currentVel);
    
    //Transform mujoco data into controller data
    state = mujocoDataToControllerData(robot, currentPos, currentVel);

    std::cout << "t = " << clock.getTime() << std::endl;
    std::cout << "Mujoco x com = " << sim.getCoM()(0) << " Controller x com = " << robot.getCoM()(0) << std::endl;
    std::cout << "Mujoco y com = " << sim.getCoM()(1) << " Controller y com = " << robot.getCoM()(1) << std::endl;
    std::cout << "Mujoco z com = " << sim.getCoM()(2) << " Controller z com = " << robot.getCoM()(2) << std::endl;  
    std::cout<< std::endl;  

    return state;
}

Eigen::MatrixXd relabelMujocoMatrix(const Robot& robot)
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

Eigen::VectorXd mujocoDataToControllerData(const Robot& robot, const Eigen::VectorXd& currentPos, const Eigen::VectorXd& currentVel)
{
    int n = robot.getNumJoints();
    Eigen::VectorXd state(2*n);

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
    
    relabelJoints = (relabelMujocoMatrix(robot).transpose()) * (currentPos.segment(7,24));
    state.segment(6,24) = relabelJoints;

    Eigen::VectorXd relabelJointsVel(24); 
    relabelJointsVel = (relabelMujocoMatrix(robot).transpose()) * (currentVel.segment(6,24));
    state.segment(n+6,24) = relabelJointsVel;

    return state;
}

void updateController(const Robot& robot, Controller& controller, Clock& clock, MujocoSim& sim)
{
    static int step_counter = 0;
    int n = robot.getNumJoints();
    Eigen::VectorXd state(2*n);

    Eigen::MatrixXd L = relabelMujocoMatrix(robot);

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

    Eigen::VectorXd tau_test = Eigen::VectorXd::Zero(24);
    Eigen::VectorXd tau0(24);
  
    state = updateState(robot,clock,sim);
      
    ControllerInput in;
    in.q    = state.segment(0,n);
    in.dq   = state.segment(n,n);
    in.time = clock.getTime();

    if (step_counter % 10 == 0) {
        controller.standStep(in);
        clock.step(); 
    }
    
    controller.inverseDynamics(in);
                  
    tau0 = controller.getTorques();

    tau_test = L*tau0; 
     
    sim.applyTorquesV2(joints, tau_test);
    step_counter++;
}
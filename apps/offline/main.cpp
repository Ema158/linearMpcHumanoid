#include "linearMpcHumanoid/robotInfo/Robot.hpp"
#include "linearMpcHumanoid/controller/controller.hpp"
#include "linearMpcHumanoid/controller/invKinematics.hpp"
#include "linearMpcHumanoid/controller/mpcLinearPendulum.hpp"
#include "linearMpcHumanoid/general/Clock.hpp"
#include "linearMpcHumanoid/general/Task.hpp"

void stand(Robot& robot, Controller& controller, Clock& clock);

Eigen::VectorXd dynamics(const Eigen::VectorXd& state, double t, Robot& Robot, Controller& controller);

int main() {
    double simulationTime = 2;
    double timeStep = 0.01;

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
    com << 0.00, 0.02, 0.26 ;
    
    //Inverse kinematics to compute the initial joint configuration
    Eigen::VectorXd desOp = ik.desiredOperationalState(nao,Rf,Lf,com);
    ik.compute(nao, desOp);
    
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
    
    //Compute the torques that balance the robot 
    stand(nao, controller, clock);

    return 0;
}

void stand(Robot& robot, Controller& controller, Clock& clock) 
{
     int n = robot.getNumJoints();
     Eigen::VectorXd state(2*n);
     state.segment(0,n) = robot.getJoints();
     state.segment(n,n) = robot.getJointsVelocity();
     
     while(std::abs(clock.getTime() - clock.getSimulationTime()) > 0.01)
     {
        //auto start = std::chrono::high_resolution_clock::now();
        
        //tau_ = out.tau;
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = end - start;
        //std::cout << "Execution time: " << elapsed.count() << " seconds\n";
        
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
        robot.setJointsVelocity(state.segment(n,n));
        clock.step();
        
     }   
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

    ControllerOutput ctrlOut = controller.standStep(in);

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
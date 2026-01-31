#pragma once
#include <Eigen/Dense>
#include <qpOASES.hpp>
//#include <chrono>

#include "linearMpcHumanoid/controller/mpcLinearPendulum.hpp"
#include "linearMpcHumanoid/controller/Dynamics.hpp"
#include "linearMpcHumanoid/controller/invKinematics.hpp"

#include "linearMpcHumanoid/robotInfo/Robot.hpp"

#include "linearMpcHumanoid/trajectories/zmpGeneration.hpp"
#include "linearMpcHumanoid/trajectories/footRefTrajectory.hpp"

#include "linearMpcHumanoid/general/Task.hpp"
#include "linearMpcHumanoid/general/Clock.hpp"
#include "linearMpcHumanoid/general/generalizedFunctions.hpp"
#include "linearMpcHumanoid/general/rk4.hpp"

/*
state = [q, v]
q = [0p1,0eta1,qJ]
    0p1 ->   base position wrt world frame
    0eta1 -> base orientation wrt world frame (Euler angles)
    qJ ->   actual joints

v = [0v1,0w1, qDJ]
    0v1 ->  linear velocity of the base (spatial velocity)
    0w1 ->  angular velocity of the base (spatial velocity)
    qDJ -> actual joints velocity
*/

struct ControllerInput {
    Eigen::VectorXd q;
    Eigen::VectorXd dq;
    double time;
};

struct ControllerOutput {
    Eigen::VectorXd tau;
};

struct WBCOutput
    {
        Eigen::VectorXd qpp;
        Eigen::VectorXd f;
        Eigen::VectorXd tau;
};

class Controller {
public:
    Controller(
        Robot& robot,
        Mpc3dLip& mpc,
        ZMP& zmp,
        std::vector<Eigen::VectorXd>& rFCoeff,
        std::vector<Eigen::VectorXd>& lFCoeff);

    ControllerOutput standStep(const ControllerInput& in);

    WBCOutput WBC(const Eigen::VectorXd& state, double t);

private: 
    Eigen::VectorXd tau_; //24 dimention vector of torques
    Eigen::VectorXd state_; // 60 dimention vector of the current configuration and velocity
    
    Robot& robot_;
    Dynamics dyn_;
    Kinematics kin_;
    Mpc3dLip& mpc_;
    ZMP zmp_;
    
    Eigen::MatrixXd frictionMatrix_ = Eigen::MatrixXd::Zero(3,4);
    double mu_ = 0.7; //friction coeff
    
    int numDynamicsEqConstraints_ = 6;
    int numFrictionEqConstraints_ = 12; //Both feet in contact with the ground
    int numFrictionIneqConstraints_ = 16 + 16; //Both feet in contact with the ground

    int numEqConstraints_ = numDynamicsEqConstraints_ + numFrictionEqConstraints_; 
    int numIneqConstraints_ = numFrictionIneqConstraints_; //ci>0, 16 for Right foot coef, 16 for left foot coef
    int numConstraints_ = numEqConstraints_ + numIneqConstraints_;
    int numDesVariables_ = 30 + 6 + 6 + 16 + 16; //joints acc(including base) + spatial force RFoot + spatial force LFoot
                                                 //...+ RFoot Coef + LFoot Coef
    int numVertex_ = 4; //Vertices in each foot
    int numCoeff_ = 4; //Number of coeff at each vertex (number of sides of pyramid friction)
    int numReactionForces_ = 12; // Number of reaction forces at current moment 12->DS 6->SS 0->noContact    
    
    //--------------------------------PD gains used in the reference accelerations-------------------------
    //PD for the joints
    Eigen::MatrixXd KpJoints_ = 300*Eigen::MatrixXd::Identity(30, 30);
    Eigen::MatrixXd KdJoints_ = 34*Eigen::MatrixXd::Identity(30, 30);

    //PD for spatial momentum
    Eigen::MatrixXd KpMom_ = 10*Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd KdMom_ = 6.32*Eigen::MatrixXd::Identity(3, 3);

    //Pd for feet position and orientation
    Eigen::MatrixXd KpFeet_ = 500*Eigen::MatrixXd::Identity(6, 6);
    Eigen::MatrixXd KdFeet_ = 44*Eigen::MatrixXd::Identity(6, 6);    
    
    const Eigen::VectorXd PDJointsAcc();
    const Eigen::VectorXd PDMomentumAcc();
    const Eigen::VectorXd PDFeetAcc(double t);

    //---------------------------------QP Weights for the WBC-----------------------------------------
    double wCoML_ = 4000; //linear momentum weight
    double wCoMK_ = 0; //angular momentum rate weight
    double wBasePos_ = 10; //base position
    double wBaseAng_ = 10; //base attitude
    double wJoints_ = 1; //rotational joints
    double wForce_ = 1; //reaction forces
    double wFoot_ = 100000; //position and orientation of both feet

    qpOASES::SQProblem qp_;
    bool qp_initialized_ = false;

    Eigen::VectorXd solveQP(qpOASES::SQProblem& qp,
    bool& initialized,
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g);

    //-----------------------------------Foot
    std::vector<Eigen::VectorXd> rFCoeff_;
    std::vector<Eigen::VectorXd> lFCoeff_;

    void frictionConstraints(Eigen::MatrixXd& Aeq,
        Eigen::VectorXd& beq,
        Eigen::MatrixXd& Aineq,
        Eigen::VectorXd& bineq);
};


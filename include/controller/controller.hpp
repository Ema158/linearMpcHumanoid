#pragma once
#include <Eigen/Dense>
#include <qpOASES.hpp>
#include "controller/Robot.hpp"
#include "controller/mpcLinearPendulum.hpp"
#include "controller/zmpGeneration.hpp"
#include "controller/Task.hpp"
#include "controller/Dynamics.hpp"
#include "controller/Clock.hpp"
#include "controller/invKinematics.hpp"
#include "controller/generalizedFunctions.hpp"
#include <memory>

/*
state = [q, v]
q = [0p,0eta,qJ]
    0p ->   base position wrt world frame
    0eta -> base orientation wrt world frame (Euler angles)
    qJ ->   actual joints

v = [0v,0w, qDJ]
    0v ->  linear velocity of the base
    0w ->  angular velocity of the base
    qDJ -> actual joints velocity
*/

class Controller {
public:
    Controller(Robot& robot, Mpc3dLip& mpc);

    void stand();

    void computeComMomentum(Eigen::VectorXd v); //Uses state and centroidal matrix to compute the velocity of the center of mass
                                //And the angular momentum of the center of mass
                               //Other option is to use the center of mass Jacobian of invKinematics

    Eigen::VectorXd WBC(double t);

    void frictionConstraints(Eigen::MatrixXd& Aeq,
        Eigen::VectorXd& beq,
        Eigen::MatrixXd& Aineq,
        Eigen::VectorXd& bineq);

private:
    double simulationTime_ = 4;
    double t_ = 0;
    double dt_ = 0.01;
    
    Eigen::VectorXd tau_; //24 dimention vector of torques
    Eigen::VectorXd state_; // 60 dimention vector of the current configuration and velocity
    
    Robot& robot_;
    Dynamics dyn_;
    Kinematics kin_;
    Mpc3dLip& mpc_;

    Eigen::Vector3d comPos_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d comVel_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d comAngMom_ = Eigen::Vector3d::Zero();
    
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
    const Eigen::VectorXd PDFeetAcc();

    //---------------------------------QP Weights for the WBC-----------------------------------------
    double wCoML_ = 4000; //linear momentum weight
    double wCoMK_ = 0; //angular momentum rate weight
    double wBasePos_ = 10; //base position
    double wBaseAng_ = 10; //base attitude
    double wJoints_ = 1; //rotational joints
    double wForce_ = 1; //reaction forces
    double wFoot_ = 100000; //position and orientation of both feet

    qpOASES::QProblem qp_;
    bool qp_initialized_ = false;

    Eigen::VectorXd solveQP(qpOASES::QProblem& qp,
    bool& initialized,
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g);
};
#include "linearMpcHumanoid/controller/controller.hpp"
#include <iostream>
#include <cmath>

Controller::Controller(Robot& robot,
    Mpc3dLip& mpc,
    ZMP& zmp,
    std::vector<Eigen::VectorXd>& rFCoeff,
    std::vector<Eigen::VectorXd>& lFCoeff
    )
    :
    robot_(robot),
    mpc_(mpc), 
    qp_(numDesVariables_,numConstraints_),
    zmp_(zmp),
    rFCoeff_(rFCoeff),
    lFCoeff_(lFCoeff)
{
    std::cout<<"Controller Initiated" << std::endl;
    std::cout<<"Initial conditions of the center of mass: " << "c = " << robot_.getCoM()(0);
    std::cout<<", "<< robot_.getCoM()(1) << ", " << robot_.getCoM()(2) << std::endl << std::endl;
    comPos_ << robot_.getCoM()(0), robot_.getCoM()(1), robot_.getCoM()(2);

    //Initial State
    state_.resize(2*robot_.getNumJoints());
    state_.segment(0,robot_.getNumJoints()) = robot_.getJoints();
    state_.segment(robot_.getNumJoints(),robot_.getNumJoints()) = robot_.getJointsVelocity();
    //Initial torques
    tau_ = Eigen::VectorXd::Zero(robot_.getNumJoints());

    //Friction matrix 
    //This matrix represent the 4 basis vectors that represent the pyramid constraint at each vertex of the foot
    frictionMatrix_.block(0,0,3,1) <<  mu_,    0, 1;
    frictionMatrix_.block(0,1,3,1) <<    0,  mu_, 1;
    frictionMatrix_.block(0,2,3,1) << -mu_,    0, 1;
    frictionMatrix_.block(0,3,3,1) <<    0, -mu_, 1;

    //QP options
    qpOASES::Options options;
    options.setToMPC();
    options.printLevel = qpOASES::PL_NONE;
    qp_.setOptions(options);

    rFCoeff_.resize(3);
    lFCoeff_.resize(3);
}

ControllerOutput Controller::standStep(const ControllerInput& in)
{
    int n = robot_.getNumJoints();

    //Update state
    robot_.updateState(in.q);
    robot_.setJointsVelocity(in.dq);

    //Update M,C, Jacobians, etc
    dyn_.computeAll(robot_);
    kin_.computeAll(robot_);

    //Update CoM state for LIP model
    Eigen::Vector2d com_xy;
    Eigen::Vector2d comVel_xy;
    com_xy << robot_.getCoM()(0), robot_.getCoM()(1);
    computeComMomentum(in.dq);
    comVel_xy << comVel_(0), comVel_(1);

    //Mpc lip model
    mpc_.compute(
        com_xy,
        comVel_xy,
        zmp_.getZmpXRef(),
        zmp_.getZmpYRef(),
        in.time
    );

    //WBC
    state_.segment(0, n) = in.q;
    state_.segment(n, n) = in.dq;
    WBCOutput out = WBC(state_, in.time);

    ControllerOutput result;
    result.tau = out.tau;
    return result;    
}

void Controller::computeComMomentum(Eigen::VectorXd v)
{
    Kinematics::swapBaseVelocityAndRefToWorldFrame(robot_.getX()[0], v);
    Eigen::VectorXd comSpatialMomentum = Eigen::VectorXd::Zero(6);

    comSpatialMomentum = dyn_.getAG()*v; 
                                                                
    comVel_ = comSpatialMomentum.segment(3,3)/robot_.getMass(); //We divided linear momentum by mass to obtain velocity 
    
    comAngMom_ = comSpatialMomentum.segment(0,3);
}

WBCOutput Controller::WBC(const Eigen::VectorXd& state, double t)
{
    int n = robot_.getNumJoints();
    Eigen::VectorXd qDD(n); // joints acceleration 
    Eigen::VectorXd tau(n-6); // torques only at actual joints, not floating base
    Eigen::VectorXd forces(numReactionForces_); // reaction forces and moments
    //Eigen::VectorXd q = Eigen::VectorXd::Zero(robot_.getNumJoints());

    //q = state.segment(0,robot_.getNumJoints());

    Eigen::VectorXd qppRef = PDJointsAcc();
    Eigen::VectorXd hGpRef = PDMomentumAcc();
    Eigen::VectorXd footAccRef = PDFeetAcc(t);

    Eigen::MatrixXd WJ = Eigen::MatrixXd::Zero(robot_.getNumJoints(), robot_.getNumJoints()); //Weight matrix of joints (including base)
    WJ.block(0,0,3,3) = wBasePos_*Eigen::Matrix3d::Identity(); // Weights base position
    WJ.block(3,3,3,3) = wBaseAng_*Eigen::Matrix3d::Identity(); // Weights base orientation
    WJ.block(6,6,robot_.getNumActualJoints(),robot_.getNumActualJoints()) 
        = wJoints_*Eigen::MatrixXd::Identity(robot_.getNumActualJoints(),robot_.getNumActualJoints()); // Weights rotational joints

    Eigen::MatrixXd WC = Eigen::MatrixXd::Zero(6,6); //Weight matrix centroidal momentum
    WC.block(0,0,3,3) = wCoMK_*Eigen::Matrix3d::Identity(); //Angular momentum
    WC.block(3,3,3,3) = wCoML_*Eigen::Matrix3d::Identity(); //Linear momentum

    Eigen::MatrixXd WFeet = wFoot_*Eigen::MatrixXd::Identity(12,12); //Weight matrix feet position and orientation

    Eigen::MatrixXd WForce = wForce_*Eigen::MatrixXd::Identity(numReactionForces_, numReactionForces_); //Weight matrix reaction spatial forces

    //Hessian construction
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(numDesVariables_, numDesVariables_); //Hessian qp
    H.block(0,0,robot_.getNumJoints(),robot_.getNumJoints()) = 
          (dyn_.getAG().transpose())*WC*(dyn_.getAG()) 
        + WJ 
        + (kin_.getFeetJacobian().transpose())*WFeet*(kin_.getFeetJacobian());

    H.block(robot_.getNumJoints(),robot_.getNumJoints(), numReactionForces_, numReactionForces_) = WForce;

    //gradient construction
    Eigen::VectorXd g = Eigen::VectorXd::Zero(numDesVariables_);
    g.segment(0,robot_.getNumJoints()) = 
              (dyn_.getAG().transpose())*WC*dyn_.getAGpqp() 
            - (dyn_.getAG().transpose())*WC*(hGpRef) 
            - WJ*qppRef
            + (kin_.getFeetJacobian().transpose())*WFeet*(dyn_.getJpqp()) 
            - (kin_.getFeetJacobian().transpose())*WFeet*(footAccRef);

    Eigen::VectorXd qpSolution = solveQP(qp_, qp_initialized_, H, g);
    forces = qpSolution.segment(n,numReactionForces_);
    Eigen::VectorXd generalizedForces(n);
    generalizedForces = (dyn_.getM())*(qpSolution.segment(0,n));// + dyn_.getC() 
                                        //- (kin_.getFeetJacobian().transpose())*x.segment(n,numReactionForces_); //tau = M*qDD + C - J*F
    tau = generalizedForces.segment(6,n-6);
    //x have the base spatial acceleration wrt to base frame, before integreated It has to be wrt world frame
    Eigen::VectorXd accBaseWrtWrld = robot_.getX()[0].colPivHouseholderQr().solve(qpSolution.segment(0,6));
    //Also the order of linear-angular base velocity is inverted in the generalized velocity vector
    qDD.segment(0,3) = accBaseWrtWrld.segment(3,3);
    qDD.segment(3,3) = accBaseWrtWrld.segment(0,3);
    qDD.segment(6,robot_.getNumActualJoints()) = qpSolution.segment(6,robot_.getNumActualJoints());

    return { qDD, forces, tau};
}

void Controller::frictionConstraints(Eigen::MatrixXd& Aeq,
    Eigen::VectorXd& beq,
    Eigen::MatrixXd& Aineq,
    Eigen::VectorXd& bineq)
{   
    //----------------------------------------Forces---------------------------------------------------------
    //Right foot
    Aeq(0,robot_.getNumJoints() + 3) = -1; //-fx + c11*v1x + c12*v2x + c13v3x + c14v4x + ... = 0
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = frictionMatrix_.row(0); 
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0); 
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0);
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(0); 

    Aeq(1,robot_.getNumJoints() + 3 + 1) = -1; //-fy + c11*v1y + c12*v2y + c13v3y + c14v4y + ... = 0
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = frictionMatrix_.row(1); 
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1); 
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1);
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(1); 

    Aeq(2,robot_.getNumJoints() + 3 + 2) = -1; //-fz + c11*v1z + c12*v2z + c13v3z + c14v4z ... = 0
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = frictionMatrix_.row(2); 
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2); 
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2);
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(2); 

    //Left foot
    Aeq(6,robot_.getNumJoints() + 3 + 6) = -1; //-fx + c11*v1x + c12*v2x + c13v3x + c14v4x + ... = 0
    Aeq.block(6,robot_.getNumJoints() + numReactionForces_ + 4*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0); 
    Aeq.block(6,robot_.getNumJoints() + numReactionForces_ + 5*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0); 
    Aeq.block(6,robot_.getNumJoints() + numReactionForces_ + 6*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0);
    Aeq.block(6,robot_.getNumJoints() + numReactionForces_ + 7*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(0); 

    Aeq(7,robot_.getNumJoints() + 3 + 6 + 1) = -1; //-fy + c11*v1y + c12*v2y + c13v3y + c14v4y + ... = 0
    Aeq.block(7,robot_.getNumJoints() + numReactionForces_ + 4*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1); 
    Aeq.block(7,robot_.getNumJoints() + numReactionForces_ + 5*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1); 
    Aeq.block(7,robot_.getNumJoints() + numReactionForces_ + 6*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1);
    Aeq.block(7,robot_.getNumJoints() + numReactionForces_ + 7*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(1); 

    Aeq(8,robot_.getNumJoints() + 3 + 6 + 2) = -1; //-fz + c11*v1z + c12*v2z + c13v3z + c14v4z ... = 0
    Aeq.block(8,robot_.getNumJoints() + numReactionForces_ + 4*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2); 
    Aeq.block(8,robot_.getNumJoints() + numReactionForces_ + 5*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2); 
    Aeq.block(8,robot_.getNumJoints() + numReactionForces_ + 6*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2);
    Aeq.block(8,robot_.getNumJoints() + numReactionForces_ + 7*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(2); 

    //----------------------------------------Moments-----------------------------------------------------------
    //Right foot
    Aeq(3,robot_.getNumJoints()) = -1; //-nx ...
    Aeq(4,robot_.getNumJoints()+1) = -1; //-ny ...
    Aeq(5,robot_.getNumJoints()+2) = -1; //-nz ...

    Eigen::MatrixXd temp;
    temp = crossMatrix(robot_.getFootVertices()[0])*frictionMatrix_; // ... + p1 x forceVertex1 + ...
     
    Aeq.block(3,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = temp.row(0);
    Aeq.block(4,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = temp.row(1);
    Aeq.block(5,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = temp.row(2);

    temp = crossMatrix(robot_.getFootVertices()[1])*frictionMatrix_; // ... + p2 x forceVertex2 + ...
    Aeq.block(3,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = temp.row(0);
    Aeq.block(4,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = temp.row(1);
    Aeq.block(5,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = temp.row(2);

    temp = crossMatrix(robot_.getFootVertices()[2])*frictionMatrix_; // ... + p3 x forceVertex3 + ...
    Aeq.block(3,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = temp.row(0);
    Aeq.block(4,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = temp.row(1);
    Aeq.block(5,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = temp.row(2);

    temp = crossMatrix(robot_.getFootVertices()[3])*frictionMatrix_; // .... + p4 x forceVertex4...
    Aeq.block(3,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1, numCoeff_) = temp.row(0);
    Aeq.block(4,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1, numCoeff_) = temp.row(1);
    Aeq.block(5,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1, numCoeff_) = temp.row(2);

    //Left foot
    Aeq(9,robot_.getNumJoints() + 6) = -1; //-nx ...
    Aeq(10,robot_.getNumJoints()+ 6 + 1) = -1; //-ny ...
    Aeq(11,robot_.getNumJoints() + 6 + 2) = -1; //-nz ...

    temp = crossMatrix(robot_.getFootVertices()[0])*frictionMatrix_; // ... + p1 x forceVertex1 + ...
     
    Aeq.block(9,robot_.getNumJoints() + numReactionForces_ + 4*numCoeff_, 1, numCoeff_) = temp.row(0);
    Aeq.block(10,robot_.getNumJoints() + numReactionForces_ + 4*numCoeff_, 1, numCoeff_) = temp.row(1);
    Aeq.block(11,robot_.getNumJoints() + numReactionForces_ + 4*numCoeff_, 1, numCoeff_) = temp.row(2);

    temp = crossMatrix(robot_.getFootVertices()[1])*frictionMatrix_; // ... + p2 x forceVertex2 + ...
    Aeq.block(9,robot_.getNumJoints() + numReactionForces_ + 5*numCoeff_, 1, numCoeff_) = temp.row(0);
    Aeq.block(10,robot_.getNumJoints() + numReactionForces_ + 5*numCoeff_, 1, numCoeff_) = temp.row(1);
    Aeq.block(11,robot_.getNumJoints() + numReactionForces_ + 5*numCoeff_, 1, numCoeff_) = temp.row(2);

    temp = crossMatrix(robot_.getFootVertices()[2])*frictionMatrix_; // ... + p3 x forceVertex3 + ...
    Aeq.block(9,robot_.getNumJoints() + numReactionForces_ + 6*numCoeff_, 1, numCoeff_) = temp.row(0);
    Aeq.block(10,robot_.getNumJoints() + numReactionForces_ + 6*numCoeff_, 1, numCoeff_) = temp.row(1);
    Aeq.block(11,robot_.getNumJoints() + numReactionForces_ + 6*numCoeff_, 1, numCoeff_) = temp.row(2);

    temp = crossMatrix(robot_.getFootVertices()[3])*frictionMatrix_; // .... + p4 x forceVertex4...
    Aeq.block(9,robot_.getNumJoints() + numReactionForces_ + 7*numCoeff_, 1, numCoeff_) = temp.row(0);
    Aeq.block(10,robot_.getNumJoints() + numReactionForces_ + 7*numCoeff_, 1, numCoeff_) = temp.row(1);
    Aeq.block(11,robot_.getNumJoints() + numReactionForces_ + 7*numCoeff_, 1, numCoeff_) = temp.row(2);

    beq = Eigen::VectorXd::Zero(12); // = 0

    //------------------------------------------Inequality-----------------------------------------
    Aineq.block(0, robot_.getNumJoints() + numReactionForces_, 2*numCoeff_*numVertex_, 2*numCoeff_*numVertex_)
                                            = Eigen::MatrixXd::Identity(2*numCoeff_*numVertex_,2*numCoeff_*numVertex_); // cij
                                            
    bineq = Eigen::VectorXd::Zero(numIneqConstraints_);                                    
}

const Eigen::VectorXd Controller::PDJointsAcc()
{
    Eigen::VectorXd qppRef = Eigen::VectorXd::Zero(robot_.getNumJoints());
    qppRef = KpJoints_*(robot_.desiredPosture() - robot_.getJoints())
     + KdJoints_*(Eigen::VectorXd::Zero(robot_.getNumJoints()) - robot_.getJointsVelocity());
    //States are x = [p, eta, qJ, v, \omega, dot{qJ}]
    //Accelerations are [\dot{\omega}, \dot{v}, \ddot{qJ}] 
    //Then we need to swap the order of the base accelerations (linear and angular)
    Eigen::VectorXd temp = qppRef.segment(0,3);
    qppRef.segment(0,3) = qppRef.segment(3,3);
    qppRef.segment(3,3) = temp;
    return qppRef;
}

const Eigen::VectorXd Controller::PDMomentumAcc()
{
    Eigen::VectorXd hGpRef = Eigen::VectorXd::Zero(6);

    Eigen::Vector3d comPosRef;
    comPosRef << mpc_.getXRef()(0), mpc_.getYRef()(0), mpc_.getZCom();
    Eigen::Vector3d comVelRef;
    comVelRef << mpc_.getXRef()(1), mpc_.getYRef()(1), 0;
    Eigen::Vector3d comAccRef;
    comAccRef << mpc_.getXRef()(2), mpc_.getYRef()(2), 0;
    
    hGpRef.segment(3,3) = robot_.getMass()*(KpMom_*(comPosRef - comPos_) + KdMom_*(comVelRef - comVel_) + comAccRef); //linear momentum
    hGpRef.segment(0,3) = KdMom_*(Eigen::Vector3d::Zero() - comAngMom_); //Zero desired angular momentum
    return hGpRef;
}

const::Eigen::VectorXd Controller::PDFeetAcc(double t)
{
    Eigen::VectorXd footAccRef = Eigen::VectorXd::Zero(12);
    Eigen::VectorXd v = robot_.getJointsVelocity();
    //Change base velocity wrt world frame -> base velocity wrt base frame
    //Beacause Jacobian is wrt base frame
    Kinematics::swapBaseVelocityAndRefToWorldFrame(robot_.getX()[0], v);
    
    //Current position of the feet
    Eigen::Vector3d rFoot = robot_.getT()[7].block(0,3,3,1);
    Eigen::Vector3d lFoot = robot_.getT()[14].block(0,3,3,1);
    
    //Compute current linear and angular velocity of the feet
    Eigen::VectorXd velR = kin_.getRightFootJacobian()*v;
    Eigen::VectorXd velL = kin_.getLeftFootJacobian()*v;

    //Desired orientation
    Eigen::Matrix3d RdesR = robot_.getRf_q0(); //Constant desired orientation parallel to ground
    Eigen::Matrix3d RdesL = robot_.getLf_q0(); //Constant desired orientation parallel to ground

    //Current Rotation matrix error in orientation
    Eigen::Matrix3d errR = RdesR.transpose()*robot_.getT()[7].block(0,0,3,3);
    Eigen::Matrix3d errL = RdesL.transpose()*robot_.getT()[14].block(0,0,3,3);

    //Current rotation axis error
    Eigen::Vector3d eR = -RdesR*rotMatToAxisAngle(errR);
    Eigen::Vector3d eL = -RdesL*rotMatToAxisAngle(errL);

    //Evaluation of polynomials
    Eigen::Vector3d rFootRef;
    Eigen::VectorXd rpFootRef(6);
    Eigen::VectorXd rppFootRef(6);

    rFootRef << polyval(rFCoeff_[0],t), polyval(rFCoeff_[1],t), polyval(rFCoeff_[2],t);
    rpFootRef << 0, 0, 0, polyval(polyder(rFCoeff_[0]),t), polyval(polyder(rFCoeff_[1]),t), polyval(polyder(rFCoeff_[2]),t);
    rppFootRef << 0, 0, 0, polyval(polyder(polyder(rFCoeff_[0])),t), polyval(polyder(polyder(rFCoeff_[1])),t), polyval(polyder(polyder(rFCoeff_[2])),t);

    Eigen::Vector3d lFootRef;
    Eigen::VectorXd lpFootRef(6);
    Eigen::VectorXd lppFootRef(6);

    lFootRef << polyval(lFCoeff_[0],t), polyval(lFCoeff_[1],t), polyval(lFCoeff_[2],t);
    lpFootRef << 0, 0, 0, polyval(polyder(lFCoeff_[0]),t), polyval(polyder(lFCoeff_[1]),t), polyval(polyder(lFCoeff_[2]),t);
    lppFootRef << 0, 0, 0, polyval(polyder(polyder(lFCoeff_[0])),t), polyval(polyder(polyder(lFCoeff_[1])),t), polyval(polyder(polyder(lFCoeff_[2])),t);

    //PD right foot
    Eigen::VectorXd rFootPosError(6);
    Eigen::VectorXd rFootVelError(6);
    rFootPosError << eR, rFootRef-rFoot;
    rFootVelError = rpFootRef - velR;
    footAccRef.segment(0,6) = KpFeet_*(rFootPosError) + KdFeet_*(rFootVelError) + rppFootRef;

    //PD left foot
    Eigen::VectorXd lFootPosError(6);
    Eigen::VectorXd lFootVelError(6);
    lFootPosError << eL, lFootRef-lFoot;
    lFootVelError = lpFootRef - velL;
    footAccRef.segment(6,6) = KpFeet_*(lFootPosError) + KdFeet_*(lFootVelError) + lppFootRef;
    return footAccRef;
}

Eigen::VectorXd Controller::solveQP(
    qpOASES::SQProblem& qp,
    bool& initialized,
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g)
{
    Eigen::VectorXd qpp(numDesVariables_);
    Eigen::MatrixXd Hsym = 0.5 * (H + H.transpose()); //guarantee symmetry

    //Dynamics constraint
    Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(numDynamicsEqConstraints_, numDesVariables_); // Dynamics constraint M /dot{v} - J^T*f = -C 
    A1.block(0,0,numDynamicsEqConstraints_,robot_.getNumJoints()) 
                    = dyn_.getM().block(0,0,numDynamicsEqConstraints_,robot_.getNumJoints());

    Eigen::MatrixXd JT = (kin_.getFeetJacobian()).transpose();
    A1.block(0,robot_.getNumJoints(),numDynamicsEqConstraints_,numReactionForces_) 
                    = -JT.block(0,0,numDynamicsEqConstraints_,numReactionForces_);
    Eigen::VectorXd b1 = -dyn_.getC().segment(0,numDynamicsEqConstraints_);

    //Friction constraints
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(numFrictionEqConstraints_,numDesVariables_); 
    Eigen::MatrixXd Aineq = Eigen::MatrixXd::Zero(numIneqConstraints_,numDesVariables_);
    Eigen::VectorXd beq = Eigen::VectorXd::Zero(numFrictionEqConstraints_); 
    Eigen::VectorXd bineq = Eigen::VectorXd::Zero(numFrictionIneqConstraints_);
    frictionConstraints(Aeq,beq,Aineq,bineq);
    
    //Express the constraints according to qpOASES documentation
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A; //CONSTRAINTS MATRIX IN QPOASES MUST BE ROW MAJOR,
                                                                            // EIGEN MATRIX IS COLUMN MAJOR BY DEFAULT
    A.setZero(numConstraints_, numDesVariables_);
    Eigen::VectorXd lbA(numConstraints_);
    Eigen::VectorXd ubA(numConstraints_);
    
    //Equality Constraints 
    //Dynamics
    A.block(0, 0, numDynamicsEqConstraints_, numDesVariables_) = A1;
    lbA.segment(0, numDynamicsEqConstraints_) = b1;
    ubA.segment(0, numDynamicsEqConstraints_) = b1;

    //Friction equality constraints
    A.block(numDynamicsEqConstraints_, 0, numFrictionEqConstraints_, numDesVariables_) = Aeq;
    lbA.segment(numDynamicsEqConstraints_, numFrictionEqConstraints_) = beq;
    ubA.segment(numDynamicsEqConstraints_, numFrictionEqConstraints_) = beq;

    //Inequality constraints 
    //Friction inequality constraints
    A.block(numDynamicsEqConstraints_ + numFrictionEqConstraints_, 0, numIneqConstraints_, numDesVariables_) = Aineq;
    lbA.segment(numDynamicsEqConstraints_ + numFrictionEqConstraints_, numIneqConstraints_) = bineq;
    ubA.segment(numDynamicsEqConstraints_ + numFrictionEqConstraints_, numIneqConstraints_).setConstant(qpOASES::INFTY); //positive infinite

    int nWSR = 500;
    qpOASES::returnValue status;
    if(!qp_initialized_){
        status = qp.init(Hsym.data(), g.data(),
            A.data(), nullptr, nullptr,
            lbA.data(), ubA.data(), nWSR);
            qp_initialized_ = true;
    }
    else{
        status = qp.hotstart(Hsym.data(), g.data(),
            A.data(), nullptr, nullptr,
            lbA.data(), ubA.data(), nWSR);
    }
    
    if (status != qpOASES::SUCCESSFUL_RETURN)
    {
    std::cerr << "QP failed, status = " << status << std::endl;
    // optionally return previous solution or zero
    }
    qp.getPrimalSolution(qpp.data());
    return qpp;
}






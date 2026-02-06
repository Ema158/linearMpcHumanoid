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
    zmp_(zmp),
    rFCoeff_(rFCoeff),
    lFCoeff_(lFCoeff), 
    qp_(numDesVariables_,numConstraints_)
    
{
    std::cout<<"Controller Initiated" << std::endl;
    std::cout<<"Initial conditions of the center of mass: " << "c = " << robot_.getCoM()(0);
    std::cout<<", "<< robot_.getCoM()(1) << ", " << robot_.getCoM()(2) << std::endl << std::endl;

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

void Controller::standStep(const ControllerInput& in)
{
    int n = robot_.getNumJoints();

    //Update state
    robot_.updateState(in.q);

    //Update M,C, Jacobians, etc
    dyn_.computeAll(robot_);
    kin_.computeAll(robot_);
    
    robot_.updateVelocityState(in.dq, dyn_.getAG());
    //Update CoM state for LIP model
    Eigen::Vector2d com_xy;
    Eigen::Vector2d comVel_xy;
    com_xy << robot_.getCoM()(0), robot_.getCoM()(1);
    comVel_xy << robot_.getComVel()(0), robot_.getComVel()(1);
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

    tau_ = out.tau;   
}

WBCOutput Controller::WBC(const Eigen::VectorXd& state, double t)
{
    int n = robot_.getNumJoints();
    Eigen::VectorXd qDD(n); // joints acceleration 
    Eigen::VectorXd torques(n-6); // torques only at actual joints, not floating base
    Eigen::VectorXd forces(numReactionForces_); // reaction forces and moments

    Eigen::VectorXd qppRef = PDJointsAcc();
    
    Eigen::VectorXd hGpRef = PDMomentumAcc();
    
    Eigen::VectorXd footAccRef = PDFeetAcc(t);
    //std::cout<<footAccRef<<std::endl;
    Eigen::MatrixXd WJ = Eigen::MatrixXd::Zero(n, n); //Weight matrix of joints (including base)
    WJ.block(0,0,3,3) = wBasePos_*Eigen::Matrix3d::Identity(); // Weights base position
    WJ.block(3,3,3,3) = wBaseAng_*Eigen::Matrix3d::Identity(); // Weights base orientation
    WJ.block(6,6,robot_.getNumActualJoints(),robot_.getNumActualJoints()) 
        = wJoints_*Eigen::MatrixXd::Identity(robot_.getNumActualJoints(),robot_.getNumActualJoints()); // Weights rotational joints

    Eigen::MatrixXd WC = Eigen::MatrixXd::Zero(6,6); //Weight matrix centroidal momentum
    WC.block(0,0,3,3) = wCoMK_*Eigen::Matrix3d::Identity(); //Angular momentum
    WC.block(3,3,3,3) = wCoML_*Eigen::Matrix3d::Identity(); //Linear momentum

    Eigen::MatrixXd WFeet = wFoot_*Eigen::MatrixXd::Identity(numReactionForces_, numReactionForces_); //Weight matrix feet position and orientation

    Eigen::MatrixXd WForce = wForce_*Eigen::MatrixXd::Identity(numReactionForces_, numReactionForces_); //Weight matrix reaction spatial forces

    //Hessian construction
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(numDesVariables_, numDesVariables_); //Hessian qp
    H.block(0,0,n,n) = 
          (dyn_.getAG().transpose())*WC*(dyn_.getAG()) 
        + WJ 
        + (kin_.getFeetJacobian().transpose())*WFeet*(kin_.getFeetJacobian());

    H.block(n, n, numReactionForces_, numReactionForces_) = WForce;

    const double eps = 1e-8;
    H.block(numDesVariablesJoints_ + numDesVariablesForces_, numDesVariablesJoints_ + numDesVariablesForces_,
        numDesVariablesCoeff_,
        numDesVariablesCoeff_)
        = eps * Eigen::MatrixXd::Identity(
        numDesVariablesCoeff_,
        numDesVariablesCoeff_);
    //gradient construction
    Eigen::VectorXd g = Eigen::VectorXd::Zero(numDesVariables_);
    
    g.segment(0, n) = 
              (dyn_.getAG().transpose())*WC*dyn_.getAGpqp() 
            - (dyn_.getAG().transpose())*WC*(hGpRef) 
            - WJ*qppRef
            + (kin_.getFeetJacobian().transpose())*WFeet*(dyn_.getJpqp()) 
            - (kin_.getFeetJacobian().transpose())*WFeet*(footAccRef);
    
    Eigen::VectorXd qpSolution = solveQP(H, g);
    
    forces = qpSolution.segment(numDesVariablesJoints_, numDesVariablesForces_);
    Eigen::VectorXd generalizedForces(n);
    generalizedForces = (dyn_.getM())*(qpSolution.segment(0,n)) + dyn_.getC() 
                                        - (kin_.getFeetJacobian().transpose())*qpSolution.segment(n,numReactionForces_); //tau = M*qDD + C - J*F
    torques = generalizedForces.segment(6,n-6);
    
    //x have the base spatial acceleration wrt to base frame, before integreated It has to be wrt world frame
    Eigen::VectorXd accBaseWrtWrld = robot_.getX()[0].colPivHouseholderQr().solve(qpSolution.segment(0,6));
    //Also the order of linear-angular base velocity is inverted in the generalized velocity vector
    qDD.segment(0,3) = accBaseWrtWrld.segment(3,3);
    qDD.segment(3,3) = accBaseWrtWrld.segment(0,3);
    qDD.segment(6,robot_.getNumActualJoints()) = qpSolution.segment(6,robot_.getNumActualJoints());

    WBCOutput out;
    out.qpp = qDD;
    out.f = forces;
    out.tau = torques;
    return out;
}

void Controller::frictionConstraints(Eigen::MatrixXd& Aeq,
    Eigen::VectorXd& beq,
    Eigen::MatrixXd& Aineq,
    Eigen::VectorXd& bineq)
{   
    Aeq.resize(numFrictionEqConstraints_, numDesVariables_);
    Aineq.resize(numFrictionIneqConstraints_, numDesVariables_);
    beq.resize(numFrictionEqConstraints_);
    bineq.resize(numFrictionIneqConstraints_);

    Aeq.setZero();
    Aineq.setZero();
    beq.setZero();
    bineq.setZero();

    const int idx_nR    = robot_.getNumJoints(); //30
    const int idx_fR    = idx_nR + 3; //33
    const int idx_nL    = idx_fR + 3; //36
    const int idx_fL    = idx_nL + 3; //39
    const int idx_muR   = idx_fL + 3; //42
    const int idx_muL   = idx_muR + numCoeff_*numVertex_; //58

    const int idx_ConstraintfR = 0; //0
    const int idx_ConstraintnR = idx_ConstraintfR + 3; //3
    const int idx_ConstraintfL = idx_ConstraintnR + 3; //6
    const int idx_ConstraintnL = idx_ConstraintfL + 3; //9
    
    //----------------------------------------Forces---------------------------------------------------------
    //Right foot
    Aeq(idx_ConstraintfR,idx_fR) = -1; //-fx + c11*v1x + c12*v2x + c13v3x + c14v4x + ... = 0
    for(int i=0;i<numVertex_;i++){
        Aeq.block(idx_ConstraintfR,idx_muR + i*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0);
    }
    
    Aeq(idx_ConstraintfR + 1,idx_fR + 1) = -1; //-fy + c11*v1y + c12*v2y + c13v3y + c14v4y + ... = 0
    for(int i=0;i<numVertex_;i++){
        Aeq.block(idx_ConstraintfR + 1,idx_muR + i*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1);
    }
    
    Aeq(idx_ConstraintfR + 2,idx_fR + 2) = -1; //-fz + c11*v1z + c12*v2z + c13v3z + c14v4z ... = 0
    for(int i=0;i<numVertex_;i++){
        Aeq.block(idx_ConstraintfR + 2,idx_muR + i*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2);
    }
    
    //----------------------------------------------------Left foot
    Aeq(idx_ConstraintfL,idx_fL) = -1; //-fx + c11*v1x + c12*v2x + c13v3x + c14v4x + ... = 0
    for(int i=0;i<numVertex_;i++){
        Aeq.block(idx_ConstraintfL, idx_muL + i*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0);
    }
    
    Aeq(idx_ConstraintfL + 1,idx_fL + 1) = -1; //-fy + c11*v1y + c12*v2y + c13v3y + c14v4y + ... = 0
    for(int i=0;i<numVertex_;i++){
        Aeq.block(idx_ConstraintfL + 1, idx_muL + i*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1);
    } 
    
    Aeq(idx_ConstraintfL + 2,idx_fL + 2) = -1; //-fz + c11*v1z + c12*v2z + c13v3z + c14v4z ... = 0
    for(int i=0;i<numVertex_;i++){
        Aeq.block(idx_ConstraintfL + 2, idx_muL + i*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2);
    }
    
    //----------------------------------------Moments-----------------------------------------------------------
    //Right foot
    Aeq(idx_ConstraintnR,idx_nR) = -1; //-nx ...

    Aeq(idx_ConstraintnR + 1,idx_nR + 1) = -1; //-ny ...

    Aeq(idx_ConstraintnR + 2, idx_nR + 2) = -1; //-nz ...

    Eigen::MatrixXd temp;
    temp = crossMatrix(robot_.getFootVertices()[0])*frictionMatrix_; // ... + p1 x forceVertex1 + ...  
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnR + i, idx_muR, 1, numCoeff_) = temp.row(i);
    }
    
    temp = crossMatrix(robot_.getFootVertices()[1])*frictionMatrix_; // ... + p2 x forceVertex2 + ...
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnR + i, idx_muR + numCoeff_, 1, numCoeff_) = temp.row(i);
    }

    temp = crossMatrix(robot_.getFootVertices()[2])*frictionMatrix_; // ... + p3 x forceVertex3 + ...
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnR + i, idx_muR + 2*numCoeff_, 1, numCoeff_) = temp.row(i);
    }

    temp = crossMatrix(robot_.getFootVertices()[3])*frictionMatrix_; // .... + p4 x forceVertex4...
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnR + i, idx_muR + 3*numCoeff_, 1, numCoeff_) = temp.row(i);
    }
    
    //Left foot
    Aeq(idx_ConstraintnL,idx_nL) = -1; //-nx ...

    Aeq(idx_ConstraintnL + 1,idx_nL + 1) = -1; //-ny ...

    Aeq(idx_ConstraintnL + 2, idx_nL + 2) = -1; //-nz ...

    temp = crossMatrix(robot_.getFootVertices()[0])*frictionMatrix_; // ... + p1 x forceVertex1 + ...
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnL + i, idx_muL, 1, numCoeff_) = temp.row(i);
    }

    temp = crossMatrix(robot_.getFootVertices()[1])*frictionMatrix_; // ... + p2 x forceVertex2 + ...
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnL + i, idx_muL + numCoeff_, 1, numCoeff_) = temp.row(i);
    }

    temp = crossMatrix(robot_.getFootVertices()[2])*frictionMatrix_; // ... + p3 x forceVertex3 + ...
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnL + i, idx_muL + 2*numCoeff_, 1, numCoeff_) = temp.row(i);
    }

    temp = crossMatrix(robot_.getFootVertices()[3])*frictionMatrix_; // .... + p4 x forceVertex4...
    for(int i=0;i<3;i++){
        Aeq.block(idx_ConstraintnL + i, idx_muL + 3*numCoeff_, 1, numCoeff_) = temp.row(i);
    }

    beq = Eigen::VectorXd::Zero(numFrictionEqConstraints_); // = 0
    //------------------------------------------Inequality-----------------------------------------
    
    Aineq.block(0, idx_muR, 2*numCoeff_*numVertex_, 2*numCoeff_*numVertex_)
                                            = Eigen::MatrixXd::Identity(2*numCoeff_*numVertex_,2*numCoeff_*numVertex_); // cij
                                            
    bineq = Eigen::VectorXd::Zero(numFrictionIneqConstraints_); 
    
    assert(Aeq.rows() == numFrictionEqConstraints_);
    assert(Aeq.cols() == numDesVariables_);
    assert(Aineq.rows() == numIneqConstraints_);
    assert(Aineq.cols() == numDesVariables_);

    if (!Aeq.allFinite()) {
    std::cerr << "Aeq has NaN or Inf!" << std::endl;
    std::abort();
    }

    if (!Aineq.allFinite()) {
        std::cerr << "Aineq has NaN or Inf!" << std::endl;
        std::abort();
    }
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
    
    hGpRef.segment(3,3) = robot_.getMass()*(KpMom_*(comPosRef - robot_.getCoM())
                         + KdMom_*(comVelRef - robot_.getComVel()) + comAccRef); //linear momentum
    hGpRef.segment(0,3) = KdMom_*(Eigen::Vector3d::Zero() - robot_.getComAngMom()); //Zero desired angular momentum
    return hGpRef;
}

const::Eigen::VectorXd Controller::PDFeetAcc(double t)
{
    Eigen::VectorXd footAccRef = Eigen::VectorXd::Zero(12);
    Eigen::VectorXd v = robot_.getJointsVelocity();
    //Change base velocity wrt world frame -> base velocity wrt base frame
    //Beacause Jacobian is wrt base frame
    swapBaseVelocityAndRefToWorldFrame(robot_.getX()[0], v);
    
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
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g)
{
    Eigen::VectorXd qpp(numDesVariables_);
    Eigen::MatrixXd Hsym = 0.5 * (H + H.transpose()); //guarantee symmetry
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Hqp;
    Hqp = Hsym;

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
    Eigen::MatrixXd Aineq = Eigen::MatrixXd::Zero(numFrictionIneqConstraints_,numDesVariables_);
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

    Hqp_ = Hsym;
    A_ = A;
    g_   = g;
    lbA_ = lbA;
    ubA_ = ubA;
    //qpOASES::SQProblem qp2(numDesVariables_, numConstraints_);

    auto check = [](const char* name, const Eigen::MatrixXd& M) {
    if (!M.allFinite()) {
        std::cerr << name << " has NaN or Inf\n";
        std::abort();
    }
    };

    auto checkv = [](const char* name, const Eigen::VectorXd& v) {
        if (!v.allFinite()) {
            std::cerr << name << " has NaN or Inf\n";
            std::abort();
        }
    };

    check("H", Hsym);
    checkv("g", g);
    check("A", A);
    checkv("lbA", lbA);
    checkv("ubA", ubA);
    status = qp_.init(Hqp_.data(), g_.data(),
        A_.data(), nullptr, nullptr,
        lbA_.data(), ubA_.data(), nWSR);
        //qp_initialized_ = true;
    
    if (status != qpOASES::SUCCESSFUL_RETURN)
    {
    std::cerr << "QP failed, status = " << status << std::endl;
    // optionally return previous solution or zero
    }
    qp_.getPrimalSolution(qpp.data()); 
    return qpp;
}






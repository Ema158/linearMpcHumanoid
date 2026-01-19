#include "controller/controller.hpp"
#include <iostream>
#include <cmath>

Controller::Controller(Robot& robot,
    Mpc3dLip& mpc)
    :
    robot_(robot),
    mpc_(mpc) 
{
    std::cout<<"Controller Initiated" << std::endl;
    std::cout<<"Initial conditions of the center of mass: " << std::endl << "c = " << robot_.getCoM();
    //Initial State
    state_.resize(2*robot_.getNumJoints());
    state_.segment(0,robot_.getNumJoints()) = robot_.getJoints();
    state_.segment(robot_.getNumJoints(),robot_.getNumJoints()) = robot_.getJointsVelocity();
    //Initial torques
    tau_ = Eigen::VectorXd::Zero(robot_.getNumJoints());

    //Friction matrix 
    //This matrix represent the 4 basis vectors that represent the pyramid constraint at each vertex of the foot
    frictionMatrix_.block(0,0,3,1) << mu_,0,1;
    frictionMatrix_.block(0,1,3,1) << 0,mu_,1;
    frictionMatrix_.block(0,2,3,1) << -mu_,0,1;
    frictionMatrix_.block(0,3,3,1) << 0,-mu_,1;

    //QP options
    qpOASES::Options options;
    options.setToMPC();
    options.printLevel = qpOASES::PL_NONE;
    qp_.setOptions(options);

}

void Controller::stand() 
{
     Clock clock;
     ZMP zmp(Task::Stand,simulationTime_,dt_,SupportFoot::Double);
     Eigen::Vector2d com_xy;
     Eigen::Vector2d comVel_xy = Eigen::Vector2d::Zero();
     Eigen::VectorXd newState = Eigen::VectorXd::Zero(6);
     com_xy << robot_.getCoM()(0),robot_.getCoM()(1);
     while(std::abs(clock.getTime() - simulationTime_) > 0.01)
     {
        dyn_.computeAll(robot_); //Computes M,C,AG,AGpqp
        kin_.computeAll(robot_); //Computes Jacobian and Jacobian bias Jpqp
        //com_xy << robot_.getCoM()(0),robot_.getCoM()(1); //Current position center of mass
        computeComMomentum(state_.segment(robot_.getNumJoints(),robot_.getNumJoints())); //Current velocity center of mass
        //comVel_xy << comVel_(0),comVel_(1);
        
        newState = mpc_.compute(com_xy,comVel_xy,zmp.getZmpXRef(),zmp.getZmpYRef(), clock.getTime());
        com_xy << newState(0),newState(3);
        comVel_xy << newState(1),newState(4);
        std::cout<<newState<<std::endl<<std::endl;

        Eigen::VectorXd qpp = WBC(clock.getTime());
        clock.step();
     }   
}

void Controller::computeComMomentum(Eigen::VectorXd v)
{
    Kinematics::swapBaseVelocityAndRefToWorldFrame(robot_.getX()[0], v);
    Eigen::VectorXd comSpatialMomentum = Eigen::VectorXd::Zero(6);

    comSpatialMomentum = dyn_.getAG()*v; 
                                                                
    comVel_ = comSpatialMomentum.segment(3,3)/robot_.getMass(); //We divided linear momentum by mass to obtain velocity 
    
    comAngMom_ = comSpatialMomentum.segment(0,3);
}

Eigen::VectorXd Controller::WBC(double t)
{
    Eigen::VectorXd qDD = Eigen::VectorXd::Zero(2*robot_.getNumJoints());
    Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(6,robot_.getNumJoints() + 12 + 16 + 16); // Dynamics constraint M /dot{v} - J^T*f = -C 
    A1.block(0,0,6,robot_.getNumJoints()) = dyn_.getM().block(0,0,6,robot_.getNumJoints());
    Eigen::MatrixXd JT = (kin_.getFeetJacobian()).transpose();
    A1.block(0,robot_.getNumJoints(),6,12) = -JT.block(0,0,6,12);
    Eigen::VectorXd b1 = -dyn_.getC().segment(0,6);

    Eigen::VectorXd qppRef = PDJointsAcc();
    Eigen::VectorXd hGpRef = PDMomentumAcc();
    Eigen::VectorXd footAccRef = PDFeetAcc();

    Eigen::MatrixXd WJ = Eigen::MatrixXd(robot_.getNumJoints(), robot_.getNumJoints()); //Weight matrix of joints (including base)
    WJ.block(0,0,3,3) = wBasePos_*Eigen::Matrix3d::Identity(); // Weights base position
    WJ.block(3,3,3,3) = wBaseAng_*Eigen::Matrix3d::Identity(); // Weights base orientation
    WJ.block(6,6,robot_.getNumActualJoints(),robot_.getNumActualJoints()) 
        = wJoints_*Eigen::MatrixXd(robot_.getNumActualJoints(),robot_.getNumActualJoints()); // Weights rotational joints

    Eigen::MatrixXd WC = Eigen::MatrixXd(6,6); //Weight matrix centroidal momentum
    WC.block(0,0,3,3) = wCoMK_*Eigen::Matrix3d::Identity(); //Angular momentum
    WC.block(3,3,3,3) = wCoML_*Eigen::Matrix3d::Identity(); //Linear momentum

    Eigen::MatrixXd WFeet = wFoot_*Eigen::MatrixXd::Identity(12,12); //Weight matrix feet position and orientation

    Eigen::MatrixXd WForce = wForce_*Eigen::MatrixXd::Identity(12,12); //Weight matrix reaction spatial forces

    //Hessian construction
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(numDesVariables_, numDesVariables_); //Hessian qp
    H.block(0,0,robot_.getNumJoints(),robot_.getNumJoints()) = 
        (dyn_.getAG().transpose())*WC*(dyn_.getAG()) + (kin_.getFeetJacobian().transpose())*WFeet*(kin_.getFeetJacobian());

    H.block(robot_.getNumJoints(),robot_.getNumJoints(),12,12) = WForce;

    //gradient construction
    Eigen::VectorXd g = Eigen::VectorXd::Zero(numDesVariables_);
    g.segment(0,robot_.getNumJoints()) = (dyn_.getAGpqp().transpose())*WC*dyn_.getAG() 
            - (hGpRef.transpose())*WC*(dyn_.getAG()) - (qppRef.transpose())*WJ 
            + (dyn_.getJpqp().transpose())*WFeet*(kin_.getFeetJacobian()) - (footAccRef.transpose())*WFeet*(kin_.getFeetJacobian());

    qDD = solveQP(qp_, qp_initialized_, H, g);
    return qDD;
}

void Controller::frictionConstraints(Eigen::MatrixXd& Aeq,
    Eigen::VectorXd& beq,
    Eigen::MatrixXd& Aineq,
    Eigen::VectorXd& bineq)
{   
    //----------------------------------------Forces---------------------------------------------------------
    Aeq(0,robot_.getNumJoints() + 3) = -1; //-fx + c11*v1x + c12*v2x + c13v3x + c14v4x + ... = 0
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = frictionMatrix_.row(0); 
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0); 
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(0);
    Aeq.block(0,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(0); 

    Aeq(0,34) = -1; //-fy + c11*v1y + c12*v2y + c13v3y + c14v4y + ... = 0
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = frictionMatrix_.row(1); 
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1); 
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(1);
    Aeq.block(1,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(1); 

    Aeq(0,35) = -1; //-fz + c11*v1z + c12*v2z + c13v3z + c14v4z ... = 0
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_, 1, numCoeff_) = frictionMatrix_.row(2); 
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_ + numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2); 
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_ + 2*numCoeff_, 1, numCoeff_) = frictionMatrix_.row(2);
    Aeq.block(2,robot_.getNumJoints() + numReactionForces_ + 3*numCoeff_, 1 ,numCoeff_) = frictionMatrix_.row(2); 

    //----------------------------------------Moments-----------------------------------------------------------
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

    beq = Eigen::VectorXd::Zero(6); // = 0

    //------------------------------------------Inequality-----------------------------------------
    Aineq.block(0, robot_.getNumJoints() + numReactionForces_, 2*numCoeff_, 2*numCoeff_)
                                            = Eigen::MatrixXd::Identity(2*numCoeff_,2*numCoeff_); // cij
                                            
    bineq = Eigen::VectorXd::Zero(2*numCoeff_);                                    
}

const Eigen::VectorXd Controller::PDJointsAcc()
{
    Eigen::VectorXd qppRef = Eigen::VectorXd::Zero(robot_.getNumJoints());
    qppRef = KpJoints_*(robot_.getJoints() - robot_.desiredPosture())
     + KdJoints_*(robot_.getJointsVelocity() - Eigen::VectorXd::Zero(robot_.getNumJoints()));
    return qppRef;
}

const Eigen::VectorXd Controller::PDMomentumAcc()
{
    Eigen::VectorXd hGpRef = Eigen::VectorXd::Zero(6);

    Eigen::Vector3d comPosRef;
    comPosRef << mpc_.getXRef()(0), mpc_.getYRef()(1), mpc_.getZCom();
    Eigen::Vector3d comVelRef;
    comVelRef << mpc_.getXRef()(1), mpc_.getYRef()(1), 0;
    Eigen::Vector3d comAccRef;
    comAccRef << mpc_.getXRef()(2), mpc_.getYRef()(2), 0;

    hGpRef.segment(3,3) = robot_.getMass()*(KpMom_*(comPos_ - comPosRef) + KdMom_*(comVel_ - comVelRef) + comAccRef); //linear momentum
    hGpRef.segment(0,3) = KdMom_*(comAngMom_ - Eigen::Vector3d::Zero()); //Zero desired angular momentum
    return hGpRef;
}

const::Eigen::VectorXd Controller::PDFeetAcc()
{
    Eigen::VectorXd footAccRef = Eigen::VectorXd::Zero(12);
    Eigen::VectorXd v = robot_.getJointsVelocity();
    //Change base velocity wrt world frame -> base velocity wrt base frame
    //Beacause Jacobian is wrt base frame
    Kinematics::swapBaseVelocityAndRefToWorldFrame(robot_.getX()[0], v);
    
    //Compute linear and angular velocity of the feet
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
    return footAccRef;
}

Eigen::VectorXd Controller::solveQP(
    qpOASES::QProblem& qp,
    bool& initialized,
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g)
{
    Eigen::MatrixXd Hsym = 0.5 * (H + H.transpose()); //guarantee symmetric and positive definite

    //Friction constraints
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(6,numDesVariables_);
    Eigen::MatrixXd Aineq = Eigen::MatrixXd::Zero(2*numCoeff_,numDesVariables_);
    Eigen::VectorXd beq = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd bineq = Eigen::VectorXd::Zero(2*numCoeff_);
    frictionConstraints(Aeq,beq,Aineq,bineq);

    //Express the constraints according to qpOASES documentation
    Eigen::MatrixXd A(numConstraints_, numDesVariables_);
    Eigen::VectorXd lbA(numConstraints_);
    Eigen::VectorXd ubA(numConstraints_);

    //Inequality constraints
    A.block(0,0,numIneqConstraints_, numDesVariables_) = Aineq;
    lbA.segment(0,numIneqConstraints_) = bineq;
    ubA.segment(0, numIneqConstraints_).setConstant(qpOASES::INFTY); //positive infinite

    //Equality Constraints 
    A.block(numIneqConstraints_,0, numEqConstraints_, numDesVariables_) = Aeq;
    lbA.segment(numIneqConstraints_, numEqConstraints_) = beq;
    ubA.segment(numIneqConstraints_, numEqConstraints_) = beq; 

    int nWSR = 50;

    if (!initialized) {
        qp.init(Hsym.data(), g.data(),
                A.data(), nullptr, nullptr,
                lbA.data(), ubA.data(), nWSR);
        initialized = true;
    } else {
        qp.hotstart(g.data(),
                    nullptr, nullptr,
                    lbA.data(), ubA.data(), nWSR);
    }

    Eigen::VectorXd qpp(numDesVariables_);
    qp.getPrimalSolution(qpp.data());
    return qpp;
}




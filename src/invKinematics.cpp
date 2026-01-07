#include "controller/invKinematics.hpp"
#include <iostream>
#define BASEDOF 6

Eigen::VectorXd invKinematics::desiredOperationalState(robotInfo robot, const Eigen::VectorXd Rf, const Eigen::VectorXd Lf, const Eigen::Vector3d com){
    Eigen::VectorXd Qd = Eigen::VectorXd::Zero(robot.getNumJoints());
    Eigen::VectorXd q = robot.getJoints();
    Qd.segment(0,6) = Rf;
    Qd.segment(6,6) = Lf;
    Qd.segment(6,12) = q.segment(12,12);
    return Qd;
}

Eigen::VectorXd invKinematics::compute(robotInfo robot, Eigen::VectorXd desOp){
    Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.getNumJoints());//Joint vector
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(robot.getNumJoints());//Actual operational variables
    Eigen::VectorXd e = Eigen::VectorXd::Zero(robot.getNumJoints());//Error vector
    Q = operationalState(robot);
    e = desOp - Q;
    double errorCriterion = e.cwiseAbs().maxCoeff(); //Absolute value of the max error in vector e
    int iter = 0; //Iteration counter
    int maxIter = 200; //Max number of iterations (solution not found)
    double tolerance = 1e-10;
    while(errorCriterion>tolerance&&iter<maxIter){ //Newton's method to solve Inverse Kinematics

    }
    return q;
}

Eigen::VectorXd invKinematics::operationalState(robotInfo robot){
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(robot.getNumJoints());
    std::vector<Eigen::Matrix4d> T = robot.getT();
    Eigen::VectorXd q = robot.getJoints(); 
    Q.segment(0,3) = T[7].block(0,3,3,1); //Position of the right foot
    //The frame attached to the foot does not coincide when the inertial frame when q=0 (rotRF(q=0) =! identity)
    //Then, if we want Roll,Pitch,Yaw angles of the foot wrt inertial frame, the rotation matrix of the foot...
    //... must change such that when q=0 rotRF(q=0) = identity
    Eigen::Matrix3d rotRF; //Rotation matrix of the right foot such that when q=0 rotRF(q=0) = identity
    rotRF = T[7].block(0,0,3,3)*robot.Rf_q0;
    Q(5) = std::atan2(rotRF(1,0),rotRF(0,0)); //Yaw Right foot
    Q(4) = std::atan2(-rotRF(2,0),std::cos(Q(5))*rotRF(0,0)+std::sin(Q(5))*rotRF(1,0)); //Pitch right foot
    Q(3) = std::atan2(std::sin(Q(5))*rotRF(0,2)-std::cos(Q(5))*rotRF(1,2),-std::sin(Q(5))*rotRF(0,1)+std::cos(Q(5))*rotRF(1,1)); //Roll right foot

    Q.segment(6,3) = T[14].block(0,3,3,1); //Position of the left foot
    Eigen::Matrix3d rotLF; //Rotation matrix of the left foot such that when q=0 rotRF(q=0) = identity
    rotLF = T[14].block(0,0,3,3)*robot.Lf_q0;
    Q(11) = std::atan2(rotLF(1,0),rotLF(0,0)); //Yaw Right foot
    Q(10) = std::atan2(-rotLF(2,0),std::cos(Q(11))*rotLF(0,0)+std::sin(Q(11))*rotLF(1,0)); //Pitch right foot
    Q(9) = std::atan2(std::sin(Q(11))*rotLF(0,2)-std::cos(Q(11))*rotLF(1,2),-std::sin(Q(11))*rotLF(0,1)+std::cos(Q(11))*rotLF(1,1));
    
    Q.segment(12,12) = q.segment(12,robot.getNumJoints()); //arms and head joints
    Q.segment(24,3) = q.segment(3,3);
    Q.segment(27,3) = robot.getCoM();
    return Q;
}

Eigen::MatrixXd invKinematics::feetJacobian(robotInfo robot){
    int dof = robot.getNumJoints(); //number of degrees of freedom
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(12,dof); //The Jacobians (6xdof) of both feet into a single Jacobian (12xdof)
    std::vector<Eigen::Matrix4d> T = robot.getT(); //Transformation matrix od each frame wrt world (0Ti)
    std::vector<Eigen::MatrixXd> X; //Velocity matrix transformation of each frame wrt to its parent
    X.resize(robot.getNumFrames());
    X = robot.getX();
    
    //We have piXi, we need nXi (8Xi for right foot and 15Xi for left foot)
    Eigen::MatrixXd JacR = frameJacobian(X,7,robot); //Frame 8 (index 7) is the sole of right foot
    Eigen::MatrixXd JacL = frameJacobian(X,14,robot); //Frame 15 (index 14) is the sole of left foot
    
    //Jacobians are wrt the local frame of the foot, they need to be rotated wrt world frame
    Eigen::MatrixXd R07 = Eigen::MatrixXd::Zero(6,6);
    R07.block(0,0,3,3) = T[7].block(0,0,3,3);
    R07.block(3,3,3,3) = T[7].block(0,0,3,3);
    JacR = R07*JacR;  

    Eigen::MatrixXd R014 = Eigen::MatrixXd::Zero(6,6);
    R014.block(0,0,3,3) = T[14].block(0,0,3,3);
    R014.block(3,3,3,3) = T[14].block(0,0,3,3);
    JacL = R014*JacL;
    //The complete Jacobian is the concatenation of the Jacobians of each foot
    J.block(0,0,6,dof) = JacR;
    J.block(6,0,6,dof) = JacL;
    
    return J;
}

Eigen::MatrixXd invKinematics::frameJacobian(std::vector<Eigen::MatrixXd> X, int frame, robotInfo robot){
    std::vector<Eigen::MatrixXd> Xn; //Velocity transformation matrix of frame n (at the foot) wrt to each frame
    std::vector<Eigen::MatrixXd> X_new; //subVector of X that only contain the kinematic chaain from base to frame n
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6,30);
    Eigen::VectorXd S = Eigen::VectorXd::Zero(6);
    S << 0,0,1,0,0,0; //Unit vector about the rotation axis of each frame (in z bc DH notation)
    std::vector<int> ant = robot.parentFrame();
    std::vector<int> act = robot.actuatedFrames();
    int numFrame = 1; //number of frames from the base to the end effector
    int i=frame;
    while(ant[i]>0){
        numFrame++;
        i--;
    }
    Xn.resize(numFrame);
    X_new.resize(numFrame);
    int j=0;
    for(int i=numFrame-1;i>=0;i--){  
        Xn[i] = Eigen::MatrixXd::Zero(6,6);
        X_new[i] = X[frame-j];
        j++;
    }
    Xn[numFrame-1] = X_new[numFrame-1];
    i = numFrame-1; 
    j = frame-1;
    do{
        J.block(0,act[j]+BASEDOF-1,6,1) = Xn[i]*S;
        Xn[ant[i]] = Xn[i]*X_new[i-1];
        i--;
        j--;
        //
    }while(ant[i+1]>0);
    J.block(0,0,6,6) = Xn[0] * Eigen::Matrix<double, 6, 6>::Identity(); //Frame 1 is a 6 dof joint
    
    return J;
}

Eigen::MatrixXd invKinematics::jacInvKinematics(robotInfo robot){
    int dof = robot.getNumJoints();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(dof,dof);
    return J;
}
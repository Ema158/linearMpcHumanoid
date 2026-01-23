#include "controller/invKinematics.hpp"
#include <iostream>
#define BASEDOF 6

Kinematics::Kinematics()
{
    JR_.resize(6,30);
    JL_.resize(6,30);
    JFeet_.resize(6,30);
}

Eigen::VectorXd Kinematics::desiredOperationalState(
    const Robot& robot,
    const Eigen::VectorXd& Rf, //Position and orientation of the right foot
    const Eigen::VectorXd& Lf, //Position and orientation of the left foot
    const Eigen::Vector3d& com)
{
    Eigen::VectorXd Qd = Eigen::VectorXd::Zero(robot.getNumJoints());
    Eigen::VectorXd q = robot.getJoints();
    Qd.segment(0,6) = Rf;
    Qd.segment(6,6) = Lf;
    Qd.segment(12,12) = q.segment(12 + BASEDOF,12);
    Qd.segment(24,3) = Eigen::Vector3d::Zero();
    Qd.segment(27,3) = com;
    return Qd;
}

void Kinematics::compute(
    Robot& robot,
    const Eigen::VectorXd& desOp)
{
    Eigen::VectorXd q = robot.getJoints();//Joint vector
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(robot.getNumJoints());//Actual operational variables
    Eigen::VectorXd e = Eigen::VectorXd::Zero(robot.getNumJoints());//Error vector
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(robot.getNumJoints(),robot.getNumJoints());
    Q = operationalState(robot);  
    e = desOp - Q;
    double errorCriterion = e.cwiseAbs().maxCoeff(); //Absolute value of the max error in vector e
    int iter = 0; //Iteration counter
    int maxIter = 200; //Max number of iterations (solution not found)
    double tolerance = 1e-10;
    while(errorCriterion>tolerance&&iter<maxIter){ //Newton's method to solve Inverse Kinematics
        J = jacInvKinematics(robot);
        q += J.colPivHouseholderQr().solve(e);
        robot.updateState(q);
        Q = operationalState(robot);
        e = desOp - Q;
        errorCriterion = e.cwiseAbs().maxCoeff();
    }
    if(iter>=maxIter){
        std::cout<<"Inv Kinematics no solution founded"<<std::endl;
    }
}

Eigen::VectorXd Kinematics::operationalState(
    const Robot& robot)
{
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(robot.getNumJoints());
    std::vector<Eigen::Matrix4d> T = robot.getT();
    Eigen::VectorXd q = robot.getJoints(); 
    Q.segment(0,3) = T[7].block(0,3,3,1); //Position of the right foot
    Q.segment(3,3) = rotMatrixToEulerAngles(T[7].block(0,0,3,3),robot.getRf_q0()); //Euler angles right foot
    
    Q.segment(6,3) = T[14].block(0,3,3,1); //Position of the left foot
    Q.segment(9,3) = rotMatrixToEulerAngles(T[14].block(0,0,3,3),robot.getLf_q0()); //Euler angles left foot
    
    Q.segment(12,12) = q.segment(12 + BASEDOF,12); //arms and head joints
    Q.segment(24,3) = q.segment(3,3);
    Q.segment(27,3) = robot.getCoM();
    return Q;
}

Eigen::MatrixXd Kinematics::feetJacobian(
    const Robot& robot)
{
    int dof = robot.getNumJoints(); //number of degrees of freedom
    Eigen::MatrixXd JFeet = Eigen::MatrixXd::Zero(12,dof); //The Jacobians (6xdof) of both feet into a single Jacobian (12xdof)
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
    JFeet.block(0,0,6,dof) = JacR;
    JFeet.block(6,0,6,dof) = JacL;

    //std::cout<<JFeet.block(0,0,6,12)<<std::endl<<std::endl;
    //std::cout<<JFeet.block(6,0,6,12)<<std::endl<<std::endl;
    return JFeet;
}

Eigen::MatrixXd Kinematics::frameJacobian(
    const std::vector<Eigen::MatrixXd>& X,
    int frame,
    const Robot& robot)
{
    std::vector<Eigen::MatrixXd> Xn; //Velocity transformation matrix of frame n (at the foot) wrt to each frame
    std::vector<Eigen::MatrixXd> X_new; //subVector of X that only contain the kinematic chain from base to frame n

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

Eigen::MatrixXd Kinematics::jacInvKinematics(
    const Robot& robot)
{
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(robot.getNumJoints(),robot.getNumJoints());
    Eigen::MatrixXd Jf = feetJacobian(robot);//Jacobian matrix of feet using spatial notation
    Eigen::VectorXd q = robot.getJoints();
    Eigen::Vector3d eta = q.segment(3,3);
    std::vector<Eigen::Matrix4d> T = robot.getT();

    //The jacobian of feetJacobian() its different from the one that we need
    //First in Jf the order of the velocity of the base is [angular velocity, linear velocity]
    //We need to change the order to [linear velocity, angular velocity] (swap columns (1:3) and (4:6))
    Eigen::MatrixXd temp = Jf.block(0,0,12,3);
    Jf.block(0,0,12,3) = Jf.block(0,3,12,3);
    Jf.block(0,3,12,3) = temp;

    //Then, Jf returns the velocity of the foot in the order [angular velocity, linear velocity]
    //We nee to change the order to [linear velocity, angular velocity] (swap rows (1:3) and (4:6))
    //right foot
    temp = Jf.block(0,0,3,robot.getNumJoints());
    Jf.block(0,0,3,robot.getNumJoints()) = Jf.block(3,0,3,robot.getNumJoints());
    Jf.block(3,0,3,robot.getNumJoints()) = temp;
    //left foot
    temp = Jf.block(6,0,3,robot.getNumJoints());
    Jf.block(6,0,3,robot.getNumJoints()) = Jf.block(9,0,3,robot.getNumJoints());
    Jf.block(9,0,3,robot.getNumJoints()) = temp;
    
    //Now Jf was computed in function of the angular velocity of the base
    //We need Jf in terms of the euler angles rate of change
    Eigen::Matrix3d Omega = matrixAngularVelToEulerDot(eta);
    //right foot
    Jf.block(0,3,3,3) = Jf.block(0,3,3,3)*Omega.inverse();
    Jf.block(3,3,3,3) = Jf.block(3,3,3,3)*Omega.inverse();
    //left foot
    Jf.block(6,3,3,3) = Jf.block(6,3,3,3)*Omega.inverse();
    Jf.block(9,3,3,3) = Jf.block(9,3,3,3)*Omega.inverse();
    
    //Finally Jf was computed to (when multiply by the joints velocities) return the angular velocity of each foot
    //We need that returns the rate of change of euler angles
    //Right foot
    Eigen::Vector3d etaFoot = rotMatrixToEulerAngles(T[7].block(0,0,3,3),robot.getRf_q0());
    Eigen::Matrix3d OmegaFoot = matrixAngularVelToEulerDot(etaFoot);
    Jf.block(3,3,3,3) = OmegaFoot*Jf.block(3,3,3,3);
    //Left foot 
    etaFoot = rotMatrixToEulerAngles(T[14].block(0,0,3,3),robot.getLf_q0());
    OmegaFoot = matrixAngularVelToEulerDot(etaFoot);
    Jf.block(9,3,3,3) = OmegaFoot*Jf.block(9,3,3,3);

    J.block(0,0,12,robot.getNumJoints()) = Jf; //Jacobian of both feet
    J.block(12,12+BASEDOF,12,12) = Eigen::MatrixXd::Identity(12,12); //Jacobian for the arms and head joints is an identity
    J.block(24,3,3,3) = Eigen::Matrix3d::Identity(3,3); //Jacobian for the base orientation is an identity
    J.block(27,0,3,robot.getNumJoints()) = comJacobian(robot); //Jacobian for the center of mass
    return J;
}

Eigen::MatrixXd Kinematics::comJacobian(
    const Robot& robot)
{
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3,robot.getNumJoints()); //center of mass jacobian whole robot
    Eigen::MatrixXd JX = Eigen::MatrixXd::Zero(3,robot.getNumJoints()); //jacobian of each center of mass  
    Eigen::Vector3d pCom = Eigen::Vector3d::Zero(3); //vector of each com wrt to world frame
    Eigen::VectorXd pComj4 = Eigen::VectorXd::Zero(4); //vector of each com wrt to local frame in homogeneous coord
    std::vector<linkInertia> links = robot.getLinks();
    std::vector<Eigen::Matrix4d> T = robot.getT();
    std::vector<Eigen::Matrix3d> crossMatrices;
    Eigen::VectorXd q = robot.getJoints();
    Eigen::Vector3d eta = q.segment(3,3);
    
    std::vector<int> ant = robot.parentFrame();
    std::vector<int> act = robot.actuatedFrames();
    crossMatrices.resize(robot.getNumFrames());
    for (int i=0;i<robot.getNumFrames()-1;i++){
        crossMatrices[i] = crossMatrix(T[i].block(0,2,3,1));

        if(links[i].mass!=0){
            pComj4 << links[i].com,1.0;
            pCom = T[i].block(0,0,3,4)*pComj4;
            JX.block(0,0,3,6) = baseJacobian(T[0].block(0,3,3,1),pCom);
            int j=i;
            while(j!=0){
                if(act[j]!=0){//only actuated frames
                    JX.block(0,act[j]+BASEDOF-1,3,1) = crossMatrices[j]*(pCom-T[j].block(0,3,3,1));
                }
                j = ant[j];
            }
            J = J + links[i].mass*JX;
            JX = Eigen::MatrixXd::Zero(3,robot.getNumJoints());
        }
    }
    J = J/robot.getMass();
    Eigen::Matrix3d Omega = matrixAngularVelToEulerDot(eta);
    J.block(0,3,3,3) = J.block(0,3,3,3)*Omega.inverse();
    return J;
}

Eigen::MatrixXd Kinematics::baseJacobian(
    const Eigen::Vector3d& pBase,
    const Eigen::Vector3d& pFrame)
{
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3,6);
    J.block(0,0,3,3) = Eigen::Matrix3d::Identity();
    J.block(0,3,3,3) = crossMatrix(pBase-pFrame);
    return J;
}

Eigen::Vector3d Kinematics::rotMatrixToEulerAngles(
    const Eigen::Matrix3d& R,
    const Eigen::Matrix3d& refFrame)
{
    Eigen::Vector3d eta = Eigen::Vector3d::Zero(3);
    Eigen::Matrix3d newR = Eigen::Matrix3d::Zero(3,3);
    newR = R*refFrame; //Change of reference of rotation matrix, such that when R=refFrame->newR=Identity->eta=0
    eta(2) = std::atan2(newR(1,0),newR(0,0)); //Yaw Right foot
    eta(1) = std::atan2(-newR(2,0),std::cos(eta(2))*newR(0,0)+std::sin(eta(2))*newR(1,0)); //Pitch right foot
    eta(0) = std::atan2(std::sin(eta(2))*newR(0,2)-std::cos(eta(2))*newR(1,2),-std::sin(eta(2))*newR(0,1)+std::cos(eta(2))*newR(1,1)); //Roll right foot
    return eta;
}

void Kinematics::swapBaseVelocityAndRefToWorldFrame(const Eigen::MatrixXd& X01, Eigen::VectorXd& v)
{
    // v is in the order [linear velocity, angular velocity]. We need to change
    // to the order [angular vel, linear velocity]
    Eigen::Vector3d temp = v.segment(0,3);
    v.segment(0,3) = v.segment(3,3);
    v.segment(3,3) = temp;

    //Also base velocity in v is wrt to world frame
    //In centroidal matrix base velocity in v is wrt base frame
    v.segment(0,6) = X01*v.segment(0,6);
}

void Kinematics::computeAll(Robot& robot)
{
    JFeet_ = feetJacobian(robot);
    JR_ = JFeet_.block(0,0,6,robot.getNumJoints());
    JL_ = JFeet_.block(6,0,6,robot.getNumJoints());
}


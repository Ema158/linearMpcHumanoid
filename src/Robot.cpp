#include "linearMpcHumanoid/robotInfo/Robot.hpp"
#define BASEDOF 6
constexpr double pi = 3.14159265358979323846;

Robot::Robot(){
    links_ = createNaoParameters();
    Eigen::VectorXd q0 = Eigen::VectorXd::Zero(getNumJoints());
    T_.resize(getNumFrames());
    updateState(q0);
    //The inertial information taken from the aldebaran documantation is wrt to an aldebaran frame at each joint
    //Aldebaran frame at each joint coincides with the orientation of the world frame when q=0
    //Algorithms suppose that inertial information is wrt joint frame
    //We need to change the inertial information to joint frame
    std::vector<Eigen::Matrix3d> Rj; //Rotation matrix of each frame
    Rj.resize(getNumFrames());
    mass_ = 0;
    for (int i=0; i< getNumFrames(); i++){
        Rj[i] = T_[i].block(0,0,3,3); //jR0
        links_[i].com = Rj[i].transpose()*links_[i].com; //jR0*0p = jp
        links_[i].inertia = Rj[i].transpose()*links_[i].inertia*Rj[i]; //jR0*0I*0Rj = jI
        mass_ += links_[i].mass;
    }
    q0 = initialConfiguration();
    updateState(q0);

    allVelocityMatrices(parentTransMatrix(T_));

    Rf_q0_ << 0,0,1, //Rotation matrix of right foot frame when q=0
             0,-1,0,
             1,0,0;
    Lf_q0_ = Rf_q0_; //Rotation matrix of left foot frame when q=0

    v_ = Eigen::VectorXd::Zero(getNumJoints()); //vector of joints velocities

    S_.resize(6);
    S_ << 0,0,1,0,0,0;

    footVertices_.resize(4);
    footVertices_[0] << 0.1, 0.025,0; //Upper right
    footVertices_[1] << 0.1, -0.025,0; //Upper left
    footVertices_[2] << -0.05, 0.025,0; //Lowe right
    footVertices_[3] << -0.05,-0.025,0; //Lower left
}

void Robot::forwardKinematics(){
    
    //std::vector<Eigen::Matrix4d> T;
    //T.resize(getNumFrames());
    Eigen::Vector3d basePos;
    basePos << q_(0), q_(1), q_(2); //Cartesian coordinates of the base of the robot
    Eigen::Vector3d baseAttitude = {q_(3), q_(4), q_(5)}; //Base attitude in Euler angles
    T_[0].block(0,3,3,1) = basePos;
    T_[0].block(0,0,3,3) = eulerAnglesToSO3(baseAttitude);
    T_[0].block(3,0,1,4) << 0,0,0,1;
    
    std::vector<double> theta;
    theta.resize(getNumActualJoints()); //theta is the Denavit Hartenberg parameter related with q_
    
    theta[0] = q_(0 + BASEDOF);
    theta[1] = q_(1 + BASEDOF) + (3.0/4)*pi;
    theta[2] = q_(2 + BASEDOF);
    theta[3] = q_(3 + BASEDOF);
    theta[4] = q_(4 + BASEDOF);
    theta[5] = q_(5 + BASEDOF);

    theta[6] = q_(6 + BASEDOF) - (1.0/2)*pi;
    theta[7] = q_(7 + BASEDOF) + (1.0/4)*pi;
    theta[8] = q_(8 + BASEDOF);
    theta[9] = q_(9 + BASEDOF);
    theta[10] = q_(10 + BASEDOF);
    theta[11] = q_(11 + BASEDOF); 
    
    theta[12] = q_(12 + BASEDOF);
    theta[13] = q_(13 + BASEDOF) + (1.0/2)*pi;
    theta[14] = q_(14 + BASEDOF);
    theta[15] = q_(15 + BASEDOF);
    theta[16] = q_(16 + BASEDOF);

    theta[17] = q_(17 + BASEDOF);
    theta[18] = q_(18 + BASEDOF) + (1.0/2)*pi;
    theta[19] = q_(19 + BASEDOF);
    theta[20] = q_(20 + BASEDOF);
    theta[21] = q_(21 + BASEDOF);

    theta[22] = q_(22 + BASEDOF);
    theta[23] = q_(23 + BASEDOF) - (1.0/2)*pi;
    theta[24] = -pi/2;
    
    std::vector<Eigen::Matrix4d> Temp = matTrans(theta); //Transformation matrix of each frame wrt its parent frame
    
    //===================Transformation matrix of the frame 1(RHipYawPitch) wrt to world frame when q_=0
    Eigen::Matrix4d auxT01;
    auxT01 <<      0, -1,      0, 0,
              0.7071,  0, 0.7071, 0,
             -0.7071,  0, 0.7071, 0,
                   0,  0,      0, 1;

    //===================Transformation matrix of the frame 9(LHipYawPitch) wrt to world frame when q_=0
    Eigen::Matrix4d auxT09;
    auxT09 << 1,       0,      0, 0,
              0,  0.7071, 0.7071, 0,
              0, -0.7071, 0.7071, 0,
              0,       0,      0, 1;

    //=================Transformation matrix of the sole of the right foot wrt to the right ankle
    Eigen::Matrix4d auxT_RFoot;
    auxT_RFoot << 1,   0,  0, -0.0452,
                  0,   1,  0,       0,
                  0,   0,  1,       0,
                  0,   0,  0,       1; 

    //=================Transformation matrix of the sole of the left foot wrt to the left ankle
    Eigen::Matrix4d auxT_LFoot;
    auxT_LFoot << 1, 0, 0, -0.0452,
                  0, 1, 0,       0,
                  0, 0, 1,       0,
                  0, 0, 0,       1; 

    //=======================Forward kinematics right leg
    T_[1] = T_[0]*auxT01*Temp[0];
    //std::cout << Temp[0] << std::endl;
    for(int i=1; i<6; i++){
        T_[i+1] = T_[i]*Temp[i]; 
    }
    T_[7] = T_[6]*auxT_RFoot;
    //======================Forward kinematics left leg
    T_[8] = T_[0]*auxT09*Temp[6];
    for(int i=8; i<13; i++){
       T_[i+1] = T_[i]*Temp[i-1]; //Temp[i-1] beacause there is now an extra frame in the sole
    }
    T_[14] = T_[13]*auxT_LFoot;

    //=====================Forward kinematics right arm
    Eigen::Vector3d rightShoulderOffset = {0,-0.098,0.13591}; //Position of the right shoulder wrt base frame
    T_[15] = Temp[12];
    T_[15].block(0,3,3,1) = T_[15].block(0,3,3,1) + rightShoulderOffset;
    T_[15] = T_[0]*T_[15];
    for(int i=15; i<19; i++){
        T_[i+1] = T_[i]*Temp[i-2];//Temp[i-2] beacause there are now two extra frames in each sole
    }

    //==================Forward kinematics left arm
    Eigen::Vector3d leftShoulderOffset = {0,0.098,0.13591}; //Position of the left shoulder wrt base frame
    T_[20] = Temp[17];
    T_[20].block(0,3,3,1) = T_[20].block(0,3,3,1) + leftShoulderOffset;
    T_[20] = T_[0]*T_[20];
    for(int i=20; i<24; i++){
        T_[i+1] = T_[i]*Temp[i-2];//Temp[i-2] beacause there are now two extra frames in each sole
    }

    //=================Forward kinematics head
    Eigen::Vector3d headOffset = {0,0,0.1615}; //Position of the head wrt base frame
    T_[25] = Temp[22];
    T_[25].block(0,3,3,1) = T_[25].block(0,3,3,1) + headOffset;
    T_[25] = T_[0]*T_[25];
    for(int i=25; i<27; i++){
      T_[i+1] = T_[i]*Temp[i-2];//Temp[i-2] beacause there are now two extra frames in each sole  
    }
        
}

const std::vector<int> Robot::parentFrame() const{
    //base is the torso of the robot
    // i frames {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28}
    std::vector<int> p_i = {-1,0,1,2,3,4,5,6,0,8,9,10,11,12,13,0,15,16,17,18,0,20,21,22,23,0,25,26};
    return p_i;
}

const std::vector<int> Robot::actuatedFrames() const{
    //base is the torso of the robot
    // i frames {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28}
    std::vector<int> a_i = {0,1,2,3,4,5,6,0,7,8,9,10,11,12,0,13,14,15,16,17,18,19,20,21,22,23,24,0};
    return a_i;
}

std::vector<Eigen::Matrix4d> matTrans(std::vector<double> theta){
    std::vector<Eigen::Matrix4d> T;
    int N = 25; // Number of frames with joints and an extra frame for the head
    T.resize(N);
    double r1 = -0.07071; //distance from the torso to the hip
    double r7 = 0.07071;//distance from the torso to the hip
    double r15 = 0.105;
    double r17 = 0.05595;
    double r20 = 0.105;
    double r22 = 0.05595;

    double d4 = -0.1; //distance from the knee to the ankle
    double d5 = -0.1029; //distance from the hip to the knee
    double d10 = -0.1; // distance from the knee to the ankle
    double d11 = -0.1029; // distance from the hip to the knee
    double d15 = -0.015;
    double d20 = -0.015;
    double d25 = 0.030;
    std::vector<double> r = {r1,0,0,0,0,0,r7,0,0,0,0,0,0,0,r15,0,r17,0,0,r20,0,r22,0,0,0};
    std::vector<double> d = {0,0,0,d4,d5,0,0,0,0,d10,d11,0,0,0,d15,0,0,0,0,d20,0,0,0,0,d25};
    std::vector<double> alpha = {0,pi/2,pi/2,0,0,-pi/2,-pi/2,-pi/2,pi/2,0,0,-pi/2,-pi/2,pi/2,pi/2,-pi/2,pi/2,pi/2,pi/2,pi/2,-pi/2,pi/2,0,-pi/2,0};

    // Transformation matrix of each frame wrt to its parent
    // Taken from W. Khalil "Modeling, Identification and Control or Robots"
    for (int i=0; i<N; i++){
        T[i](0,0) =  std::cos(theta[i]);
        T[i](0,1) = -std::sin(theta[i]);
        T[i](0,2) = 0;
        T[i](0,3) = d[i];

        T[i](1,0) = std::cos(alpha[i])*std::sin(theta[i]);
        T[i](1,1) = std::cos(alpha[i])*std::cos(theta[i]);
        T[i](1,2) = -std::sin(alpha[i]);
        T[i](1,3) = -r[i]*std::sin(alpha[i]);

        T[i](2,0) = std::sin(alpha[i])*std::sin(theta[i]);
        T[i](2,1) = std::sin(alpha[i])*std::cos(theta[i]);
        T[i](2,2) = std::cos(alpha[i]);
        T[i](2,3) = r[i]*std::cos(alpha[i]);

        T[i](3,0) = 0;
        T[i](3,1) = 0;
        T[i](3,2) = 0;
        T[i](3,3) = 1;
    }
    
    return T;
}

void Robot::computeCoM(){
    Eigen::Vector3d com = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d pComj; //Position of the center of mass of the jth body
    Eigen::Vector4d pComj4 = Eigen::VectorXd::Zero(4); //Position of the center of mass of the jth body in homogenous coordinates 
                                                     // a 4th dimention vector with a one at the end [com,1]
    
    for (int i=0;i<getNumFrames(); i++){
        pComj4 << links_[i].com,1.0; //com of each body wrt to body frame
        
        pComj = T_[i].block(0,0,3,4)*pComj4; //com of each body wrt world frame
        com = com + links_[i].mass*pComj; //com of full robot wrt world frame
    }
    CoM_ = com/mass_;
}



Eigen::VectorXd initialConfiguration(){
    Eigen::VectorXd q = Eigen::VectorXd::Zero(30); //Initial configuration of the robot
    q << -0.0185, 0, 0.282, 0, 0, 0, //base position and orientation
        0, 0, -0.5, 0.8, -0.3, 0, // right leg
        0, 0, -0.5, 0.8, -0.3, 0, //left leg
        1.6, 0, 0, 0, 0, //right arm
        -1.6, 0, 0, 0, 0, //left arm
        0,0; // head
    return q;
}

const Eigen::VectorXd Robot::desiredPosture(){
    Eigen::VectorXd qDes = Eigen::VectorXd::Zero(30);
    qDes << -0.0185,0.00,0.282,0,0,0, //base
            0,0,-0.5,0.8,-0.3,0, //right leg
            0,0,-0.5,0.8,-0.3,0, //left leg
            1.6,0,0,0,0, //right arm
            -1.6,0,0,0,0, //left arm
            0,0; //head
    return qDes;
}

void Robot::updateState(const Eigen::VectorXd& q_new){
    q_ = q_new;
    forwardKinematics();
    computeCoM();
    allVelocityMatrices(parentTransMatrix(T_));
}

std::vector<Eigen::Matrix4d> Robot::parentTransMatrix(
    const std::vector<Eigen::Matrix4d>& T)
{
    std::vector<Eigen::Matrix4d> piTi;
    piTi.resize(getNumFrames());
    std::vector<int> ant = parentFrame();
    piTi[0] = T[0]; //Not used 
    for (int i=1;i<getNumFrames();i++){
        piTi[i] = inverseTransformationMatrix(T[ant[i]])*T[i];
    }
    return piTi;
}

void Robot::allVelocityMatrices(
    const std::vector<Eigen::Matrix4d>& piTi)
{
    std::vector<int> act = actuatedFrames();
    X_.resize(getNumFrames());
    X_[0] = velocityMatrix(piTi[0]); //Velocity matrix of the base frame
    for(int i=0;i<getNumFrames();i++){
            X_[i] = velocityMatrix(piTi[i]);
    }
}

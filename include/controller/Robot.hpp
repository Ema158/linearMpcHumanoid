// Creates the robot
#pragma once
#include <Eigen/Dense>
#include <cmath>
#include "linkInertia.hpp"
#include "robotParameters.hpp"
#include "generalizedFunctions.hpp"
#define NUM_JOINTS 30 //24 real joints plus 6 of the floating base
#define NUM_ACTUAL_JOINTS 24 // 24 rotational joints
#define NUM_FRAMES 28 //1 for base, 24 for each joint, 2 extra feet, 1 extra head
#define NUM_BODIES 25 // 

class Robot{
public:
    Robot();
    const std::vector<int> parentFrame() const; // vector that contains the parent frame p(i) of each frame i
    
    const std::vector<int> actuatedFrames() const; //vector that contains a zero if the frame is not actuated and the number of the joint if it is actuated
    
    int getNumFrames() const {return numFrames;}
    
    int getNumJoints() const {return numJoints;}
    
    int getNumActualJoints() const {return numActualJoints;}
    
    int getNumBodies() const {return numBodies;}
    
    const Eigen::VectorXd& getJoints() const {return q_;}
    
    const Eigen::VectorXd& getJointsVelocity() const {return v_;}
    
    const std::vector<Eigen::Matrix4d>& getT() const {return T_;} 
    
    const std::vector<linkInertia>& getLinks() const {return links_;}
    
    const Eigen::Vector3d& getCoM() const {return CoM_;}
    
    double getMass() const {return mass_;}
    
    Eigen::Matrix3d getRf_q0() const {return Rf_q0_;}
    
    Eigen::Matrix3d getLf_q0() const {return Lf_q0_;}
    
    const std::vector<Eigen::MatrixXd>& getX() const {return X_;}

    Eigen::VectorXd getS() const {return S_;}
    
    void computeCoM();
    
    void forwardKinematics();
    
    void updateState(const Eigen::VectorXd& q_new);

    void setJointsVelocity(const Eigen::VectorXd& v_new) {v_ = v_new;}

    std::vector<Eigen::Matrix4d> parentTransMatrix(const std::vector<Eigen::Matrix4d>& T);

    void allVelocityMatrices(const std::vector<Eigen::Matrix4d>& piTi);

    const std::vector<Eigen::Vector3d>& getFootVertices() const {return footVertices_;} 
 
    const Eigen::VectorXd desiredPosture();
private:
    int numJoints = NUM_JOINTS;
    int numActualJoints = NUM_ACTUAL_JOINTS;
    int numFrames = NUM_FRAMES;
    int numBodies = NUM_BODIES;
    Eigen::VectorXd q_; //Generalized coordinates
    Eigen::VectorXd v_; //Generalized velocities (note v \neq dot{q})
    std::vector<Eigen::Matrix4d> T_; //Transformation matrix of each frame wrt world frame
    std::vector<linkInertia> links_;
    std::vector<Eigen::MatrixXd> X_; //Velocity matrix of each frame wrt to its parent
    Eigen::Vector3d CoM_;
    double mass_;
    Eigen::Matrix3d Rf_q0_; //Rotation matrix of right foot frame when q=0
    Eigen::Matrix3d Lf_q0_; //Rotation matrix of left foot frame when q=0 
    Eigen::VectorXd S_; //screw axis joint rotation
    std::vector<Eigen::Vector3d> footVertices_;  
};

Eigen::VectorXd initialConfiguration();

std::vector<Eigen::Matrix4d> matTrans(std::vector<double> theta);


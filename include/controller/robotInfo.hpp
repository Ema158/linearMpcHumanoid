// Creates the robot
#pragma once
#include <Eigen/Dense>
#include <cmath>
#include "linkInertia.hpp"
#include "robotParameters.hpp"
#define NUM_JOINTS 30 //24 real joints plus 6 of the floating base
#define NUM_ACTUAL_JOINTS 24 // 24 rotational joints
#define NUM_FRAMES 28 //1 for base, 24 for each joint, 2 extra feet, 1 extra head
#define NUM_BODIES 25 // 

class robotInfo{
public:
    robotInfo();
    std::vector<int> parentFrame(); // vector that contains the parent frame p(i) of each frame i
    std::vector<int> actuatedFrames(); //vector that contains a zero if the frame is not actuated and the number of the joint if it is actuated
    
    int getNumFrames() const {return numFrames;}
    int getNumJoints() const {return numJoints;}
    int getNumActualJoints() const {return numActualJoints;}
    int getNumBodies() const {return numBodies;}
    Eigen::VectorXd getJoints() const {return q;}
    Eigen::VectorXd getJointsVel() const {return v;}
    std::vector<Eigen::Matrix4d> getT() const {return T;} 
    std::vector<linkInertia> getLinks() const {return links;}
    
    Eigen::VectorXd initialConfiguration();
    Eigen::Matrix3d eulerAnglesToSO3(const Eigen::Vector3d& eulerAngles);
    std::vector<Eigen::Matrix4d> forwardKinematics(Eigen::VectorXd q);
    Eigen::VectorXd desiredPosture();
    
    void setJoints(const Eigen::VectorXd q_new){q = q_new;}
    void setVelocities(const Eigen::VectorXd v_new){v = v_new;}
    void setT(const std::vector<Eigen::Matrix4d> T_new){T = T_new;}
    void setLinks(const std::vector<linkInertia> links_new){links = links_new;}

    Eigen::Matrix3d Rf_q0; //Rotation matrix of right foot frame when q=0
    Eigen::Matrix3d Lf_q0; //Rotation matrix of left foot frame when q=0 

    double mass;
    Eigen::Vector3d getCoM();
    
private:
    int numJoints = NUM_JOINTS;
    int numActualJoints = NUM_ACTUAL_JOINTS;
    int numFrames = NUM_FRAMES;
    int numBodies = NUM_BODIES;
    Eigen::VectorXd q; //Generalized coordinates
    Eigen::VectorXd v; //Generalized velocities (note v \neq dot{q})
    std::vector<Eigen::Matrix4d> T;
    std::vector<Eigen::Matrix4d> matTrans(std::vector<double> theta);
    std::vector<linkInertia> links;

};

#pragma once
#include <Eigen/Dense>
#include "controller/Task.hpp"

class ZMP{
public:
    ZMP(const Task task); //Regulation DS

    ZMP(
        const Task task,
        const double simulationTime,
        const double timeStep,
        const SupportFoot supportFoot); //Regulation DS or SS
    
    ZMP(
        const Task task,
        const int numSteps,
        const double timePerStep,
        const double simulationTime); //Walk
    
    void stanceZMP(); //ZMP reference for a standing motion 
    void walkZMP(); //ZMP reference for a walking motion

    const Eigen::VectorXd getZmpXRef() const {return zmpXRef_;}
    const Eigen::VectorXd getZmpYRef() const {return zmpYRef_;}
private:
    Eigen::VectorXd zmpXRef_;
    Eigen::VectorXd zmpYRef_;
    Eigen::VectorXd zmpXMax_; //Constraints are not neccesary since the WBC guaranteees ZMP inside support polygon
    Eigen::VectorXd zmpYMax_; //Constraints are not neccesary since the WBC guaranteees ZMP inside support polygon
    Eigen::VectorXd zmpXMin_; //Constraints are not neccesary since the WBC guaranteees ZMP inside support polygon
    Eigen::VectorXd zmpYMin_; //Constraints are not neccesary since the WBC guaranteees ZMP inside support polygon
    double simulationTime_ = 1;
    double timeStep_ = 0.01;
    SupportFoot supportFoot_ = SupportFoot::Double;
    Task task_ = Task::Stand;
    int numSteps_ = 1;
    double timePerStep_ = 0.5;
};
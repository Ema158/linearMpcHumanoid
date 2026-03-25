#pragma once
#include <Eigen/Dense>
#include "linearMpcHumanoid/general/Task.hpp"
#include "linearMpcHumanoid/controller/ContactState.hpp"

class ZMP{
public:

    ZMP(const double timeStep,
        const Contact& supportFoot); 
    
    void stanceZMP(const double simulationTime); //ZMP reference for a standing motion 
    void walkZMP(const GaitParameters& gaitParameters); //ZMP reference for a walking motion
    void updateWalkZMP();

    const Eigen::VectorXd getZmpXRef() const {return zmpXRef_;}
    const Eigen::VectorXd getZmpYRef() const {return zmpYRef_;}

    void setContact(const Contact& newContact){contact_ = newContact;};
    void setGaitParameters(const GaitParameters& newGaitParameters){gaitParameters_ = newGaitParameters;};
private:
    Eigen::VectorXd zmpXRef_;
    Eigen::VectorXd zmpYRef_;

    double simulationTime_ = 1;
    double timeStep_ = 0.01;

    Contact contact_ = Contact::Both;

    int samplesPerStep_;
    int samplesDS_;
    int samples_;

    GaitParameters gaitParameters_;
};
#include "linearMpcHumanoid/trajectories/zmpGeneration.hpp"
#include <iostream>

ZMP::ZMP(
    const double timeStep,
    const Contact& contact) 
    :
    timeStep_(timeStep),
    contact_(contact)
{
     
}

void ZMP::stanceZMP(double simulationTime)
{
    samples_ = static_cast<int>((simulationTime + 0.5) / timeStep_);
    zmpXRef_.resize(samples_);
    zmpYRef_.resize(samples_);

    zmpXRef_.setZero();
    
    switch (contact_) {
        case Contact::Right:
            zmpYRef_.setConstant(-0.05);
            break;
        case Contact::Left:
            zmpYRef_.setConstant(0.05);
            break;
        case Contact::Both:
            zmpYRef_.setZero();
            break;
    }

}

void ZMP::walkZMP(const GaitParameters& gaitParameters)
{
    samplesPerStep_ = static_cast<int>(gaitParameters.timePerStep / timeStep_);
    samplesDS_ = static_cast<int>(0.2*samplesPerStep_);
    samples_ = 2*samplesPerStep_; //A reference of two steps is generated

    zmpXRef_.resize(samples_);
    zmpYRef_.resize(samples_);

    zmpXRef_.setZero();
    zmpYRef_.setZero();

    updateWalkZMP();
}

void ZMP::updateWalkZMP()
{
    if(contact_ == Contact::Both){
        zmpYRef_.segment(0,samplesDS_) = Eigen::VectorXd::LinSpaced(samplesDS_, gaitParameters_.currentYPos, gaitParameters_.futureYPos);
        zmpYRef_.segment(samplesDS_,samplesPerStep_) = Eigen::VectorXd::Constant(samplesPerStep_, gaitParameters_.futureYPos);
        zmpYRef_.segment(samplesPerStep_,samplesPerStep_) = Eigen::VectorXd::Constant(samplesPerStep_, gaitParameters_.futureYPos);
    }

    else if(contact_ == Contact::Right || contact_ == Contact::Left){
        zmpYRef_.segment(0,samplesPerStep_) = Eigen::VectorXd::Constant(samplesPerStep_, gaitParameters_.currentYPos);
        zmpYRef_.segment(samplesPerStep_,samplesPerStep_) = Eigen::VectorXd::Constant(samplesPerStep_, gaitParameters_.futureYPos);
        
        zmpXRef_.segment(0,samplesPerStep_) = Eigen::VectorXd::Constant(samplesPerStep_, gaitParameters_.currentXPos);
        zmpXRef_.segment(samplesPerStep_,samplesPerStep_) = Eigen::VectorXd::Constant(samplesPerStep_, gaitParameters_.futureXPos);

    }   
    
}

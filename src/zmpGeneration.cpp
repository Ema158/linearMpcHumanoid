#include "linearMpcHumanoid/trajectories/zmpGeneration.hpp"
#include <iostream>

ZMP::ZMP(
    const Task task)
    :
    task_(task)
{
    stanceZMP();
}

ZMP::ZMP(
    const Task task,
    const double simulationTime,
    const double timeStep,
    const SupportFoot supportFoot) 
    :
    simulationTime_(simulationTime),
    timeStep_(timeStep),
    supportFoot_(supportFoot)
{
    stanceZMP();
}

ZMP::ZMP(
    const Task task,
    const int numSteps,
    const double timePerStep,
    const double simulationTime)
    :
    task_(task),
    numSteps_(numSteps),
    timePerStep_(timePerStep),
    simulationTime_(simulationTime)
{

}

void ZMP::stanceZMP()
{
    int samples = static_cast<int>((simulationTime_ + 0.5) / timeStep_);
    
    zmpXRef_.resize(samples);
    zmpYRef_.resize(samples);

    zmpXRef_.setZero();
    
    //std::cout<< zmpXRef_ <<std::endl;
    switch (supportFoot_) {
        case SupportFoot::Right:
            zmpYRef_.setConstant(-0.05);
            break;
        case SupportFoot::Left:
            zmpYRef_.setConstant(0.05);
            break;
        default:
            zmpYRef_.setZero();
            break;
    }

}

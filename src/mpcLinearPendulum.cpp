#include "controller/mpcLinearPendulum.hpp"
#include <iostream>

Mpc3dLip::Mpc3dLip()
{
    initialize();
}

Mpc3dLip::Mpc3dLip(
    const double timeStep,
    const double timeHorizon,
    const double zCom)
    : timeStep_(timeStep),
      timeHorizon_(timeHorizon),
      zCom_(zCom)
{
    initialize();
}

void Mpc3dLip::initialize()
{
    int horizon = static_cast<int>(timeHorizon_/timeStep_);
    Eigen::Matrix2d A_power = Eigen::Matrix2d::Identity();
    A_ << 1,timeStep_,
            0,1;
    B_ << (timeStep_*timeStep_)/2, timeStep_;
    C_ << 1,0;
    D_ = zCom_/gravity_;

    Px_ = Eigen::MatrixXd::Zero(horizon+1,2);
    Pu_ = Eigen::MatrixXd::Zero(horizon+1,horizon+1);

    Px_.block(0,0,1,2) = C_;
    Pu_(0,0) = D_;

    for(int i=1; i<=horizon;i++){
        A_power *= A_;
        Px_.block(i,0,1,2) = C_*A_power;
        Pu_(i,i-1) = C_*B_;
        Pu_(i,i) = D_;

        Eigen::Matrix2d A_j = Eigen::Matrix2d::Identity();
        for (int j=1; j<=horizon-i; j++){ //fill each column of Pu
            A_j *= A_;
            Pu_(i+j,i-1) = C_*A_j*B_;
        }
    }
}
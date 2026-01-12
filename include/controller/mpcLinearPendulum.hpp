#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

class Mpc3dLip{

public:
    Mpc3dLip();
    
    Mpc3dLip(
        const double timeStep,
        const double timeHorizon,
        const double zCoM);

    Eigen::Vector3d getXRef() const {return xRef_;}
    Eigen::Vector3d getYRef() const {return yRef_;}

    
private:
    Eigen::Vector3d xRef_; //Pos,Vel,Acc of the com in x
    Eigen::Vector3d yRef_; //Pos,Vel,Acc of the com in y
    Eigen::MatrixXd Px_;
    Eigen::MatrixXd Pu_;
    Eigen::Matrix2d A_;
    Eigen::Vector2d B_;
    Eigen::RowVector2d C_;
    double D_;    
    double timeStep_ = 0.01;
    double timeHorizon_ = 0.5;
    double zCom_ = 0.26;
    double gravity_ = 9.81;

    void initialize();
};
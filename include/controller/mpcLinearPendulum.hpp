#pragma once
#include <Eigen/Dense>
#include <qpOASES.hpp>

class Mpc3dLip{

public:
    Mpc3dLip();
    
    Mpc3dLip(
        const double dt,
        const double timeHorizon,
        const double zCoM);

    Mpc3dLip(
        const double dt,
        const double timeHorizon,
        const double zCoM,
        const double alpha,
        const double beta);

    Eigen::Vector3d getXRef() const {return xRef_;}
    Eigen::Vector3d getYRef() const {return yRef_;}

    Eigen::VectorXd compute(
        const Eigen::Vector2d& posCom,
        const Eigen::Vector2d& velCom,
        const Eigen::VectorXd& zmpXRef,
        const Eigen::VectorXd& zmpYRef,
        double t);

    
private:
    Eigen::Vector3d xRef_; //Pos,Vel,Acc of the com in x
    Eigen::Vector3d yRef_; //Pos,Vel,Acc of the com in y
    Eigen::MatrixXd Px_;
    Eigen::MatrixXd Pu_;
    Eigen::Matrix2d A_;
    Eigen::Vector2d B_;
    Eigen::RowVector2d C_;
    double D_;    
    double dt_ = 0.01;
    double timeHorizon_ = 0.5;
    double zCom_ = 0.26;
    double gravity_ = 9.81;
    double alpha_ = 1e-3;
    double beta_ = 1;
    int horizon_ = static_cast<int>(timeHorizon_ / dt_);
    qpOASES::QProblem qpx_;
    bool qpx_initialized_ = false;
    qpOASES::QProblem qpy_;
    bool qpy_initialized_ = false;

    void initialize();
};

double solve1DQP(
    qpOASES::QProblem& qp,
    bool& initialized,
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g);
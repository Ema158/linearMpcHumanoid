#include "linearMpcHumanoid/controller/mpcLinearPendulum.hpp"

Mpc3dLip::Mpc3dLip()
    : qpx_(horizon_ + 1, 0),
      qpy_(horizon_ + 1, 0)
{
    initialize();
}

Mpc3dLip::Mpc3dLip(
    const double dt,
    const double timeHorizon,
    const double zCom)
    : dt_(dt),
      timeHorizon_(timeHorizon),
      zCom_(zCom),
      qpx_(horizon_ + 1, 0),
      qpy_(horizon_ + 1, 0)

{
    initialize();
}

Mpc3dLip::Mpc3dLip(
    const double dt,
    const double timeHorizon,
    const double zCom,
    const double alpha,
    const double beta)
    : dt_(dt),
      timeHorizon_(timeHorizon),
      zCom_(zCom),
      alpha_(alpha),
      beta_(beta),
      qpx_(horizon_ + 1, 0),
      qpy_(horizon_ + 1, 0)
{
    initialize();
}

void Mpc3dLip::initialize()
{
    horizon_ = static_cast<int>(timeHorizon_/dt_);
    Eigen::Matrix2d A_power = Eigen::Matrix2d::Identity();
    A_ << 1,dt_,
            0,1;
    B_ << (dt_*dt_)/2, dt_;
    C_ << 1,0;
    D_ = -zCom_/gravity_;

    Px_ = Eigen::MatrixXd::Zero(horizon_+1,2);
    Pu_ = Eigen::MatrixXd::Zero(horizon_+1,horizon_+1);

    Px_.block(0,0,1,2) = C_;
    Pu_(0,0) = D_;

    for(int i=1; i<=horizon_;i++){
        A_power *= A_;
        Px_.block(i,0,1,2) = C_*A_power;
        Pu_(i,i-1) = C_*B_;
        Pu_(i,i) = D_;

        Eigen::Matrix2d A_j = Eigen::Matrix2d::Identity();
        for (int j=1; j<=horizon_-i; j++){ //fill each column of Pu
            A_j *= A_;
            Pu_(i+j,i-1) = C_*A_j*B_;
        }
    }

    //QP options
    qpOASES::Options options;
    options.setToMPC();
    options.printLevel = qpOASES::PL_NONE;
    qpx_.setOptions(options);
    qpy_.setOptions(options);
}

void Mpc3dLip::compute(
    const Eigen::Vector2d& posCom,
    const Eigen::Vector2d& velCom,
    const Eigen::VectorXd& zmpXRef,
    const Eigen::VectorXd& zmpYRef,
    double t)
{
    Eigen::VectorXd Xref = Eigen::VectorXd(6); //pos,vel and acc of the com in x-y
    Eigen::Vector2d xk(posCom(0), velCom(0));
    Eigen::Vector2d yk(posCom(1), velCom(1));

    Eigen::MatrixXd H = alpha_*Eigen::MatrixXd::Identity(horizon_ + 1,horizon_ + 1)
                      + beta_*(Pu_.transpose()*Pu_);

    int k = static_cast<int>(t/ dt_);
    auto zmpX_horizon = zmpXRef.segment(k, horizon_ + 1);
    auto zmpY_horizon = zmpYRef.segment(k, horizon_ + 1);

    Eigen::VectorXd gx = beta_*Pu_.transpose()*(Px_*xk - zmpX_horizon);
    Eigen::VectorXd gy = beta_*Pu_.transpose()*(Px_*yk - zmpY_horizon);

    // === QP solve ===
    double accx = solve1DQP(qpx_, qpx_initialized_, H, gx);
    double accy = solve1DQP(qpy_, qpy_initialized_, H, gy);
    
    //Compute the center of mass position and velocity
    xk = A_*xk + B_*accx;
    yk = A_*yk + B_*accy;
    
    xRef_ << xk,accx;
    yRef_ << yk,accy;
}

double solve1DQP(
    qpOASES::QProblem& qp,
    bool& initialized,
    const Eigen::MatrixXd& H,
    const Eigen::VectorXd& g)
{
    int nWSR = 50;

    if (!initialized) {
        qp.init(H.data(), g.data(),
                nullptr, nullptr, nullptr,
                nullptr, nullptr, nWSR);
        initialized = true;
    } else {
        qp.hotstart(g.data(),
                    nullptr, nullptr,
                    nullptr, nullptr, nWSR);
    }

    Eigen::VectorXd u(H.rows());
    qp.getPrimalSolution(u.data());
    return u(0);
}
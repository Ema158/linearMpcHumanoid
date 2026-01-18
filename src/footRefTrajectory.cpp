#include "controller/footRefTrajectory.hpp"

std::vector<Eigen::VectorXd> footCoeffTrajectory(Eigen::VectorXd currentPos,
    Eigen::VectorXd desPos,
    double stepHeight,
    double T)
{
    std::vector<Eigen::VectorXd> Coeff;
    Coeff.resize(3); // x, y and z trajectories

    Eigen::MatrixXd Posx;
    Posx << 0, currentPos(0),
            T, desPos(0);
    Eigen::MatrixXd Velx;
    Velx << 0, 0,
            T, 0;
    Eigen::MatrixXd Accx;
    Accx << 0, 0,
            T, 0;

    Eigen::MatrixXd Posy;
    Posy << 0, currentPos(1),
            T, desPos(1);
    Eigen::MatrixXd Vely;
    Vely = Velx;
    Eigen::MatrixXd Accy;
    Accy = Accx;

    Eigen::MatrixXd Posz;
    Posz << 0, currentPos(2),
            T/2, stepHeight;
            T, desPos(2);
    Eigen::MatrixXd Velz;
    Velz << 0, 0,
            T/2, 0,
            T, 0;
    Eigen::MatrixXd Accz;
    Accz = Accx;

    return Coeff;
}
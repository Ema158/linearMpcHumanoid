#include "linearMpcHumanoid/trajectories/footRefTrajectory.hpp"
#include <iostream>

std::vector<Eigen::VectorXd> footCoeffTrajectory(const Eigen::Vector3d& currentPos,
    const Eigen::Vector3d& desPos,
    double stepHeight,
    double T)
{
    std::vector<Eigen::VectorXd> Coeff;
    Coeff.resize(3); // x, y and z trajectories
    
    Eigen::MatrixXd Posx(2,2);
    Posx << 0, currentPos(0),
            T, desPos(0);
    
    Eigen::MatrixXd Velx(2,2);
    Velx << 0, 0,
            T, 0;
    Eigen::MatrixXd Accx(2,2);
    Accx << 0, 0,
            T, 0;
    
    Coeff[0] = findPolyCoeff(Posx,Velx,Accx);
    
    Eigen::MatrixXd Posy(2,2);
    Posy << 0, currentPos(1),
            T, desPos(1);
    Eigen::MatrixXd Vely(2,2);
    Vely = Velx;
    Eigen::MatrixXd Accy(2,2);
    Accy = Accx;
    Coeff[1] = findPolyCoeff(Posy,Vely,Accy);
        
    Eigen::MatrixXd Posz(3,2);
    Posz << 0, currentPos(2),
            T/2, stepHeight,
            T, desPos(2);
    Eigen::MatrixXd Velz(3,2);
    Velz << 0, 0,
            T/2, 0,
            T, 0;
    Eigen::MatrixXd Accz(2,2);
    Accz = Accx;
    Coeff[2] = findPolyCoeff(Posz,Velz,Accz);
        
    return Coeff;
}
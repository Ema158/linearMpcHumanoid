#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/Robot.hpp"
#include "controller/invKinematics.hpp"
#include "controller/dynamics.hpp"
#include <iostream>

int main() {
    Robot nao;
    Eigen::VectorXd Rf = Eigen::VectorXd::Zero(6);
    Rf(1) = -0.05;
    Eigen::VectorXd Lf = Eigen::VectorXd::Zero(6);
    Lf(1) = 0.05;
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    com(2) = 0.26;
    Eigen::VectorXd desOp = ik::desiredOperationalState(nao,Rf,Lf,com);
    //std::cout<< desOp << std::endl;
    Eigen::VectorXd q = ik::compute(nao, desOp);
    std::cout<<q<<std::endl;
    //invKinematics invK;
    //dynamics dyn;
    
    /*Eigen::MatrixXd J = Eigen::MatrixXd::Zero(12,30);
    J = ik::jacInvKinematics(nao);

    std::cout<<J.block(0,0,6,12) << std::endl << std::endl;
    std::cout<<J.block(0,12,6,12) << std::endl << std::endl;
    std::cout<<J.block(6,0,6,12) << std::endl << std::endl;
    std::cout<<J.block(6,12,6,12) << std::endl << std::endl;*/

    //std::vector<Eigen::MatrixXd> I = dyn.allSpatialInertiaMatrices(nao);
    //Eigen::VectorXd C = dyn.computeC(nao,I,true);
    //std::cout<<C<<std::endl;
    //Eigen::MatrixXd M = dyn.computeM(nao,I);
    //Eigen::MatrixXd AG = dyn.centroidalMatrix(M,nao);

    //std::cout<<AG.block(0,0,6,12)<<std::endl<<std::endl;
    //std::cout<<AG.block(0,12,6,12)<<std::endl<<std::endl;
    //std::cout<<AG.block(0,24,6,6)<<std::endl<<std::endl;
    /*Eigen::MatrixXd Jcom = invK.comJacobian(nao);
    std::cout<<Jcom.block(0,0,3,12)<<std::endl<<std::endl;
    std::cout<<Jcom.block(0,12,3,12)<<std::endl<<std::endl;
    std::cout<<Jcom.block(0,24,3,6)<<std::endl<<std::endl;*/

    /*for(int i=0;i<nao.getNumFrames()-1;i++){
        std::cout<<i<<std::endl<<T[i]<<std::endl<<std::endl;
    }*/
    return 0;
}
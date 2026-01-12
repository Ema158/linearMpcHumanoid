#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/Robot.hpp"
#include "controller/invKinematics.hpp"
#include "controller/dynamics.hpp"
#include <iostream>

int main() {
    Robot nao;
    invKinematics invK;
    //dynamics dyn;
    
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(12,30);
    J = invK.jacInvKinematics(nao);

    std::cout<<J.block(0,0,6,12) << std::endl << std::endl;
    std::cout<<J.block(0,12,6,12) << std::endl << std::endl;
    std::cout<<J.block(6,0,6,12) << std::endl << std::endl;
    std::cout<<J.block(6,12,6,12) << std::endl << std::endl;

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
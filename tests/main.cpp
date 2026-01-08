#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/robotInfo.hpp"
#include "controller/invKinematics.hpp"
#include "controller/dynamics.hpp"
#include <iostream>

int main() {
    robotInfo nao;
    invKinematics invK;
    dynamics dyn;
    Eigen::VectorXd q = Eigen::VectorXd::Zero(nao.getNumJoints());
    //Desired initial position
    Eigen::VectorXd desRf = Eigen::VectorXd::Zero(6);//Desired position for the right foot
    desRf << 0,-0.05,0,0,0,0;
    Eigen::VectorXd desLf = Eigen::VectorXd::Zero(6);//Desired position and orientation for the left foot
    desLf << 0,0.05,0,0,0,0;
    Eigen::Vector3d desCom;//Desires position for the com
    desCom << 0,0,0.26;
    std::vector<Eigen::Matrix4d> T = nao.getT();
    q = nao.getJoints();
    
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(12,30);
    J = invK.feetJacobian(nao);

    /*std::cout<<J.block(0,0,6,12) << std::endl << std::endl;
    std::cout<<J.block(0,12,6,12) << std::endl << std::endl;
    std::cout<<J.block(6,0,6,12) << std::endl << std::endl;
    std::cout<<J.block(6,12,6,12) << std::endl << std::endl;*/

    std::vector<Eigen::MatrixXd> I = dyn.allSpatialInertiaMatrices(nao);
    //Eigen::VectorXd C = dyn.computeC(nao,I,true);
    //std::cout<<C<<std::endl;
    Eigen::MatrixXd M = dyn.computeM(nao,I);
    return 0;
}
#include <iostream>
#include "controller/controller.hpp"
#include "controller/robotParameters.hpp"
#include "controller/robotInfo.hpp"
#include "controller/invKinematics.hpp"
#include <iostream>

int main() {
    robotInfo nao;
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
    std::cout<< nao.getCoM() << std::endl;
    /*for (int i=0; i<nao.getNumFrames();i++){
        std::cout<< T[i] << std::endl;
    }*/
    return 0;
}
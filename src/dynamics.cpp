#include <controller/dynamics.hpp>
#include <iostream>
#define BASEDOF 6

Eigen::MatrixXd dynamics::spatialInertiaMatrix(linkInertia link){
    Eigen::MatrixXd I = Eigen::MatrixXd::Zero(6,6);
    I.block(0,0,3,3) = link.inertia - link.mass*crossMatrix(link.com)*crossMatrix(link.com);
    I.block(0,3,3,3) = link.mass*crossMatrix(link.com);
    I.block(3,0,3,3) = -link.mass*crossMatrix(link.com);
    I.block(3,3,3,3) = link.mass*Eigen::Matrix3d::Identity();
    return I;
}

std::vector<Eigen::MatrixXd> dynamics::allSpatialInertiaMatrices(robotInfo robot){
    std::vector<linkInertia> links = robot.getLinks();
    std::vector<Eigen::MatrixXd> I;
    std::vector<int> act = robot.actuatedFrames();
    I.resize(robot.getNumFrames());
    I[0] = spatialInertiaMatrix(links[0]);
    //std::cout<< I[0] << std::endl << std::endl;
    for(int i=1;i<robot.getNumFrames();i++){
        if (act[i]!=0){ //Only compute the spatial inertia of frames with actual joints
            I[i] = spatialInertiaMatrix(links[i]);
        }
    }
    return I;
}

Eigen::VectorXd dynamics::computeC(robotInfo robot, std::vector<Eigen::MatrixXd> I, bool isGravity){
   Eigen::VectorXd C = Eigen::VectorXd::Zero(robot.getNumJoints());
   Eigen::VectorXd g = Eigen::VectorXd::Zero(6); //gravity acc of the base wrt to world frame
   g(5) = isGravity*9.81;
   
   std::vector<Eigen::MatrixXd> X = robot.getX();
   Eigen::VectorXd qD = robot.getJointsVelocity();
   std::vector<int> act = robot.actuatedFrames();
   std::vector<int> ant = robot.parentFrame();
   
   Eigen::VectorXd S = Eigen::VectorXd::Zero(6);
   S << 0,0,1,0,0,0;

   std::vector<Eigen::VectorXd> vel; //spatial velocity of each frame
   vel.resize(robot.getNumFrames());
   std::vector<Eigen::VectorXd> acc; //spatial acceleration of each frame
   acc.resize(robot.getNumFrames());
   std::vector<Eigen::VectorXd> f; //spatial force of each frame
   f.resize(robot.getNumFrames());

   //Base velocity and spatial force
   vel[0] = qD.segment(0,6);
   acc[0] = X[0]*g; // transform the acc wrt world frame to wrt base frame 
   f[0] = I[0]*acc[0] + spatialCrossMatrixForce(vel[0])*I[0]*vel[0];
   //Forward pass Newton-Euler
   //std::cout<<f[0]<<std::endl<<std::endl;    
   for (int i=1;i<robot.getNumFrames();i++){
        if (act[i]!=0){
            vel[i] = X[i]*vel[ant[i]] + S*qD[act[i] + BASEDOF - 1];
            acc[i] = X[i]*acc[ant[i]] + spatialCrossMatrix(vel[i])*S*qD[act[i] + BASEDOF - 1];
            f[i] = I[i]*acc[i] + spatialCrossMatrixForce(vel[i])*I[i]*vel[i];
            //std::cout<<i<<std::endl<<f[i]<<std::endl<<std::endl;
        }
    }
    //Backward pass
    //std::cout<<std::endl<<f[26]<<std::endl<<std::endl;
    for(int i=robot.getNumFrames()-1;i>0;i--){
        //std::cout<<act[i]<<std::endl;
        if (act[i]!=0){
            C(act[i] + BASEDOF - 1) = (S.transpose())*f[i];
            std::cout<<C(act[i] + BASEDOF - 1)<<std::endl;
            f[ant[i]] = f[ant[i]] + X[i].transpose()*f[i];
            
        }
    }
   return C; 
}
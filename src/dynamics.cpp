#include <controller/dynamics.hpp>
#include <iostream>
#define BASEDOF 6

Eigen::MatrixXd dynamics::spatialInertiaMatrix(const linkInertia& link){
    Eigen::MatrixXd I = Eigen::MatrixXd::Zero(6,6);
    I.block(0,0,3,3) = link.inertia - link.mass*crossMatrix(link.com)*crossMatrix(link.com);
    I.block(0,3,3,3) = link.mass*crossMatrix(link.com);
    I.block(3,0,3,3) = -link.mass*crossMatrix(link.com);
    I.block(3,3,3,3) = link.mass*Eigen::Matrix3d::Identity();
    return I;
}

std::vector<Eigen::MatrixXd> dynamics::allSpatialInertiaMatrices(const Robot& robot){
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

Eigen::VectorXd dynamics::computeC(const Robot& robot, const std::vector<Eigen::MatrixXd>& I, bool isGravity){
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
   for (int i=1;i<robot.getNumFrames();i++){
        if (act[i]!=0){
            vel[i] = X[i]*vel[ant[i]] + S*qD[act[i] + BASEDOF - 1];
            acc[i] = X[i]*acc[ant[i]] + spatialCrossMatrix(vel[i])*S*qD[act[i] + BASEDOF - 1];
            f[i] = I[i]*acc[i] + spatialCrossMatrixForce(vel[i])*I[i]*vel[i];
        }
    }
    //Backward pass
    for(int i=robot.getNumFrames()-1;i>0;i--){
        if (act[i]!=0){
            C(act[i] + BASEDOF - 1) = (S.transpose())*f[i];
            f[ant[i]] = f[ant[i]] + X[i].transpose()*f[i];          
        }
    }
    C.segment(0,6) = f[0]; //Base joint is a 6x6 dof, thats why for base joint S = identity
   return C; 
}

Eigen::MatrixXd dynamics::computeM(const Robot& robot, const std::vector<Eigen::MatrixXd>& I){
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(robot.getNumJoints(),robot.getNumJoints()); //inertia Matrix floating base
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(robot.getNumActualJoints(),robot.getNumActualJoints()); //inertia Matrix actuated joints
    std::vector<Eigen::MatrixXd> X = robot.getX(); // 
    std::vector<int> act = robot.actuatedFrames();
    std::vector<int> ant = robot.parentFrame();
    Eigen::VectorXd S = Eigen::VectorXd::Zero(6); //axis joint 
    S << 0,0,1,0,0,0; //It is assumed all joints are rotational joints that rotates around z axis
    std::vector<Eigen::MatrixXd> Ic = I; //Ic is the vector of composed spatial inertias, 
                                         //It is initialized as the spatial inertia of each body
    std::vector<Eigen::VectorXd> f; //spatial force at each actuated frame
    f.resize(robot.getNumFrames());
    Eigen::MatrixXd F2 = Eigen::MatrixXd::Zero(6,robot.getNumActualJoints());
    //F2 is a block matrix of M
    //F2 is the effect of each unit acc joint projected in the base frame 
    int j=0;
    for (int i=robot.getNumFrames()-1;i>=0;i--){
        if (act[i]!=0){ //only considered actuated frames
            Ic[ant[i]] = Ic[ant[i]] + X[i].transpose()*Ic[i]*X[i]; //computation of composite spatial inertias
            f[i] = Ic[i]*S;//spatial force at current frame
            H(act[i]-1,act[i]-1) = S.transpose()*f[i]; //torque of the joint in the current frame
            //Now we need to project the spatial force at current frame into each frame of the current kimematic chain
            j=i;
            while(ant[j]!=0){ //while we dont reach the base
                f[i] = X[j].transpose()*f[i]; //project the spatial force into each j frame
                j = ant[j]; //go one frame back
                H(act[j]-1,act[i]-1) = S.transpose()*f[i];// torque at frame j
                H(act[i]-1,act[j]-1) = H(act[j]-1,act[i]-1); //symmetry in M matrix
            }
            F2.block(0,act[i]-1,6,1) = X[j].transpose()*f[i]; //project spatial force at frame i into base frame
                                                        //here S=identity beacause base frame is a 6 dof virtual joint
        }  
    }
    M.block(0,0,6,6) = Ic[0];
    M.block(BASEDOF,BASEDOF,robot.getNumActualJoints(),robot.getNumActualJoints()) = H;
    M.block(0,BASEDOF,6,robot.getNumActualJoints()) = F2;
    M.block(BASEDOF,0,robot.getNumActualJoints(),6) = F2.transpose();
    return M;
}

Eigen::MatrixXd dynamics::centroidalMatrix(const Eigen::MatrixXd& M, const Robot& robot){
    Eigen::MatrixXd AG = Eigen::MatrixXd(6,robot.getNumJoints()); //6 for linear and angular momentum of the com
    Eigen::MatrixXd T01 = robot.getT()[0];
    Eigen::MatrixXd Ic1 = M.block(0,0,6,6); //composed spatial inertia matrix of the base
    Eigen::MatrixXd F = M.block(0,6,6,robot.getNumActualJoints());
    Eigen::MatrixXd X1G = Eigen::MatrixXd::Zero(6,6); //velocity transformation matrix base frame/com frame
    Eigen::Vector3d p1G = Eigen::Vector3d::Zero(3);
    p1G << Ic1(2,4)/robot.getMass(),Ic1(0,5)/robot.getMass(),Ic1(1,3)/robot.getMass();
    X1G.block(0,0,3,3) = T01.block(0,0,3,3); //Rotation matrix of base frame wrt world frame
    X1G.block(3,3,3,3) = T01.block(0,0,3,3); //Rotation matrix of base frame wrt world frame
    X1G.block(0,3,3,3) = -T01.block(0,0,3,3)*crossMatrix(p1G);

    AG.block(0,0,6,6) = X1G*Ic1;
    AG.block(0,6,6,robot.getNumActualJoints()) = X1G*F;
    return AG;                                                 
}
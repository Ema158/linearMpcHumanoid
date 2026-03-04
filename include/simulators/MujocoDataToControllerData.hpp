#pragma once

#include "simulators/MujocoSim.hpp"

#include "linearMpcHumanoid/robotInfo/Robot.hpp"
#include "linearMpcHumanoid/general/Clock.hpp"
#include "linearMpcHumanoid/controller/controller.hpp"

#include <iostream>
#include <Eigen/Dense>

Eigen::MatrixXd relabelMujocoMatrix(const Robot& robot);

Eigen::MatrixXd fromMujocoFrameToControllerFrame(const Eigen::MatrixXd& R01,const Eigen::Vector3d& p01);

Eigen::VectorXd mujocoDataToControllerData(const Robot& robot, const Eigen::VectorXd& currentPos, const Eigen::VectorXd& currentVel);

Eigen::VectorXd updateState(const Robot& robot, Clock& clock, MujocoSim& sim);

void updateController(const Robot& robot, Controller& controller, Clock& clock, MujocoSim& sim);


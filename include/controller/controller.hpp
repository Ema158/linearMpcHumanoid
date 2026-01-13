#pragma once
#include <Eigen/Dense>
#include "controller/Robot.hpp"
#include "controller/mpcLinearPendulum.hpp"
#include "controller/zmpGeneration.hpp"
#include "controller/Task.hpp"

class Controller {
public:
    Controller(Robot& robot);

    void stand();

private:
    double simulationTime_ = 4;
    double t_ = 0;
    double dt_ = 0.01;
    Eigen::VectorXd tau_;
    Robot robot_;
};
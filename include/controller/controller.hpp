#pragma once
#include <Eigen/Dense>

class Controller {
public:
    Controller();

    void setState(float c,
                  float d);

    void compute();

    const int torques() const;

};
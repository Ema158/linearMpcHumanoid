#pragma once

#include <mujoco/mujoco.h>
#include <GLFW/glfw3.h>

#include "simulators/MujocoSim.hpp"

#include "linearMpcHumanoid/general/Clock.hpp"

class MujocoViewer {
public:
    MujocoViewer(MujocoSim& sim);
    ~MujocoViewer();

    void run(std::function<void()> control_cb, Clock& clock);

private:
    MujocoSim& sim_;

    GLFWwindow* window_;

    mjvScene scn_;
    mjvCamera cam_;
    mjvOption opt_;
    mjrContext con_;

    void initGLFW();
    void initMujoco();
};
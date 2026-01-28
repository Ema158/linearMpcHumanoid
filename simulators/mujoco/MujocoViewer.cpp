#include "simulators/MujocoViewer.hpp"

#include <stdexcept>
#include <iostream>

MujocoViewer::MujocoViewer(MujocoSim& sim)
    : sim_(sim), window_(nullptr)
{
    initGLFW();
    initMujoco();
}

MujocoViewer::~MujocoViewer()
{
    mjr_freeContext(&con_);
    mjv_freeScene(&scn_);
    glfwDestroyWindow(window_);
    glfwTerminate();
}

void MujocoViewer::initGLFW()
{
    if (!glfwInit())
        throw std::runtime_error("GLFW init failed");

    window_ = glfwCreateWindow(1200, 900, "MuJoCo Viewer", nullptr, nullptr);
    if (!window_)
        throw std::runtime_error("Failed to create GLFW window");

    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1); // vsync
}

void MujocoViewer::initMujoco()
{
    mjv_defaultScene(&scn_);
    mjv_makeScene(sim_.model(), &scn_, 2000);

    mjv_defaultCamera(&cam_);
    cam_.distance = 1.5;
    cam_.elevation = -20;
    cam_.azimuth = 90;

    mjv_defaultOption(&opt_);

    mjr_defaultContext(&con_);
    mjr_makeContext(sim_.model(), &con_, mjFONTSCALE_150);
}

void MujocoViewer::run(std::function<void()> control_cb, Clock& clock)
{
    
    while (!glfwWindowShouldClose(window_)&&std::abs(clock.getTime() - clock.getSimulationTime()) > 0.01) {
        if (control_cb) {
            control_cb();     //  controller called here
        }

        // Step simulation
        sim_.step();

        // Update scene
        mjv_updateScene(
            sim_.model(),
            sim_.data(),
            &opt_,
            nullptr,
            &cam_,
            mjCAT_ALL,
            &scn_
        );

        // Get framebuffer size
        int width, height;
        glfwGetFramebufferSize(window_, &width, &height);
        mjrRect viewport = {0, 0, width, height};

        // Render
        mjr_render(viewport, &scn_, &con_);

        glfwSwapBuffers(window_);
        glfwPollEvents();
    }
}

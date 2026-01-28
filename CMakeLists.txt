cmake_minimum_required(VERSION 3.16)

project(whole_body_controller LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

find_package(glfw3 REQUIRED)

# qpOASES
include_directories($ENV{HOME}/.local/include)
link_directories($ENV{HOME}/.local/lib)

# -------------------------------
# MuJoCo (prebuilt shared library)
# -------------------------------
add_library(mujoco SHARED IMPORTED GLOBAL)

set_target_properties(mujoco PROPERTIES
    IMPORTED_LOCATION "$ENV{HOME}/.mujoco/mujoco/build/lib/libmujoco.so"
    IMPORTED_NO_SONAME TRUE
    INTERFACE_INCLUDE_DIRECTORIES "$ENV{HOME}/.mujoco/mujoco/include"
)

add_library(controller_core
    src/controller.cpp
    src/robotParameters.cpp
    src/Robot.cpp
    src/invKinematics.cpp
    src/generalizedFunctions.cpp
    src/Dynamics.cpp
    src/mpcLinearPendulum.cpp
    src/zmpGeneration.cpp
    src/footRefTrajectory.cpp
    src/Clock.cpp
)

target_include_directories(controller_core PUBLIC
    include
)

target_link_libraries(controller_core PUBLIC
    Eigen3::Eigen
    qpOASES
)

# -------------------------------
# Simulators
# -------------------------------
add_library(mujoco_sim
    simulators/mujoco/MujocoSim.cpp
    simulators/mujoco/MujocoViewer.cpp
)

target_include_directories(mujoco_sim PUBLIC
    include
)

target_link_libraries(mujoco_sim PUBLIC
    controller_core
    mujoco
    glfw
)

# -------------------------------
# Executables
# -------------------------------
add_executable(offline_controller
    apps/offline/main.cpp
)
target_link_libraries(offline_controller controller_core)

add_executable(mujoco_controller
    apps/mujoco/main.cpp
)
target_link_libraries(mujoco_controller mujoco_sim)
cmake_minimum_required(VERSION 3.16)

project(whole_body_controller LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

# qpOASES
include_directories($ENV{HOME}/.local/include)
link_directories($ENV{HOME}/.local/lib)

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
)

target_include_directories(controller_core PUBLIC
    include
)

target_link_libraries(controller_core PUBLIC
    Eigen3::Eigen
    qpOASES
)

add_executable(offline_controller
    simulators/offline/main.cpp
)

target_link_libraries(offline_controller
    controller_core
)
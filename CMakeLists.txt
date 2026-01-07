cmake_minimum_required(VERSION 3.16)

project(whole_body_controller LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

add_library(controller_core
    src/controller.cpp
    src/robotParameters.cpp
    src/robotInfo.cpp
    src/invKinematics.cpp
    src/generalizedFunctions.cpp
    src/dynamics.cpp
)

target_include_directories(controller_core PUBLIC
    include
)

target_link_libraries(controller_core PUBLIC
    Eigen3::Eigen
)

add_executable(test_controller
    tests/main.cpp
)

target_link_libraries(test_controller
    controller_core
)
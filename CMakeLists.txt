cmake_minimum_required(VERSION 3.16)

project(whole_body_controller LANGUAGES CXX)

function(enable_strict_warnings target)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE
            -O1
            -g
            -Wall
            -Wextra
            -Wpedantic
            -Wshadow
            -Wconversion
            -Wsign-conversion
        )
    endif()
endfunction()

function(enable_sanitizers target)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
        )
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined
        )
    endif()
endfunction()

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

enable_strict_warnings(controller_core)
enable_sanitizers(controller_core)

target_include_directories(controller_core PUBLIC
    include
)

target_link_libraries(controller_core PUBLIC
    Eigen3::Eigen
    qpOASES
)

get_target_property(EIGEN_INC Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)

target_include_directories(controller_core SYSTEM PUBLIC
    ${EIGEN_INC}
)

target_compile_definitions(controller_core PRIVATE
    EIGEN_INITIALIZE_MATRICES_BY_NAN
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

enable_sanitizers(mujoco_sim)

# -------------------------------
# Executables
# -------------------------------
add_executable(offline_controller
    apps/offline/main.cpp
)
target_link_libraries(offline_controller controller_core)
enable_sanitizers(offline_controller)

add_executable(mujoco_controller
    apps/mujoco/main.cpp
)
target_link_libraries(mujoco_controller mujoco_sim)

enable_strict_warnings(mujoco_controller)
enable_sanitizers(mujoco_controller)
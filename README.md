# readme

⚠️ Work in progress

This repository contains an experimental linear MPC / WBC framework
for humanoid robots. The current offline simulator integrates the
closed-loop system using RK4. MuJoCo integration is planned.

This project has as objective the C++ implementation of a model based Whole-Body Controller on a floating-base humanoid robot using a linear MPC formulation and centroidal dynamics. Outputs of the controller are torques that are directly send to the robot. The results are validated in Webots in tasks such as balance, walking and jumping.

The project aims to use the fewest possible number of libraries in order to deeply understand what is happening inside each component. The only component that is not fully developed here is the qp solver.

The controller is divided in three parts: A generator of desired trajectories, the WBC that generates joint accelerations and reaction forces that tracks these trajectories such that constraints inherent to balance are fullfiled, and an inverse dynamics block that transforms these accelerations into torques that are send to the robot.

This controller was first developed on MATLAB and the reference codes are included in the MATLAB folder.

## Desired trajectories generator 
The key variables of the robot that we need to impose desired trajectories are: the position of the center of mass, the angular momentum around the center of mass, the position and orientation of the base of the robot (the torso of the robot), the position and orientation of each foot. 

Almost every desired trajectory for each key variable will be arbitrary using a polynomial trajectory in time. The CoM trajectory will be found using a linear inverted pendulum model with mpc over an horizon.

Desired trajectories are pass trough a PD controller in order to generate better reference to track by the Whole Body Controller.

## Whole Body Control
As Whole-Body Controller a qp is formulated with the following objective: minimize the error of:
 1) The rate of change of the spatial momentum (linear and angular) of the center of mass. 
 2) The acceleration of each foot.
 3) The acceleration of each joint.

The Whole-Body Control has as outputs the acceleration of the joints (including the floating base) and the reaction forces. 

## Inverse Dynamics
Using the dynamic model of the floating base robot, joints accelerations and reaction forces are converted into torques that are send to the robot.

## State estimation
An state estimator is pending to implement. Current simulations use information given by the simulator (Webots) to feedback the state of the base of the robot at each time step.




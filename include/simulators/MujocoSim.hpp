#pragma once

#include <string>
#include <Eigen/Dense>

// MuJoCo
#include <mujoco/mujoco.h>

// Your project types
#include "linearMpcHumanoid/controller/controller.hpp"   // ControllerInput, ControllerOutput

class MujocoSim {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // ---- construction / destruction ----
  MujocoSim(const std::string& model_path, const Eigen::VectorXd& q0);
  ~MujocoSim();

  // non-copyable
  MujocoSim(const MujocoSim&) = delete;
  MujocoSim& operator=(const MujocoSim&) = delete;

  // ---- simulation control ----
  void step();                 // one MuJoCo step
  void reset();                // reset to initial state

  // ---- controller interface ----
  ControllerInput getControllerInput() const;
  void applyTorques(const Eigen::VectorXd& tau);
  void applyJointPositions(const Eigen::VectorXd& qdes);

  // ---- gravity compensation ----
  Eigen::VectorXd computeGravityTorques();

  // ---- accessors ----
  int nq() const { return nq_; }
  int nv() const { return nv_; }
  int nu() const { return nu_; }

  const mjModel* model() const { return m_; }
  mjData* data() const { return d_; }

  void getMujocoState(Eigen::VectorXd& q, Eigen::VectorXd& v);

private:
  // ---- MuJoCo core ----
  mjModel* m_ = nullptr;
  mjData*  d_ = nullptr;

  // ---- model dimensions (cached) ----
  int nq_ = 0;   // number of generalized positions
  int nv_ = 0;   // number of generalized velocities
  int nu_ = 0;   // number of actuators

  // ---- helpers ----
  void loadModel(const std::string& model_path, const Eigen::VectorXd& q0);

  //-----PD gains for position controller----------------
  Eigen::VectorXd kp_, kd_;
  
};
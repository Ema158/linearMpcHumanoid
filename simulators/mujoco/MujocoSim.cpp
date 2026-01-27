#include "simulators/MujocoSim.hpp"

#include <iostream>
#include <stdexcept>

// -----------------------------
// Constructor / Destructor
// -----------------------------

MujocoSim::MujocoSim(const std::string& model_path) {
  loadModel(model_path);
}

MujocoSim::~MujocoSim() {
  if (d_) {
    mj_deleteData(d_);
    d_ = nullptr;
  }
  if (m_) {
    mj_deleteModel(m_);
    m_ = nullptr;
  }
}

// -----------------------------
// Model loading
// -----------------------------

void MujocoSim::loadModel(const std::string& model_path) {
  char error[1024] = {0};

  m_ = mj_loadXML(model_path.c_str(), nullptr, error, sizeof(error));
  if (!m_) {
    throw std::runtime_error(
      std::string("MuJoCo model load failed: ") + error
    );
  }

  d_ = mj_makeData(m_);
  if (!d_) {
    mj_deleteModel(m_);
    m_ = nullptr;
    throw std::runtime_error("MuJoCo mj_makeData failed");
  }

  // cache dimensions
  nq_ = m_->nq;
  nv_ = m_->nv;
  nu_ = m_->nu;

  std::cout << "[MujocoSim] Model loaded\n"
            << "  nq = " << nq_ << "\n"
            << "  nv = " << nv_ << "\n"
            << "  nu = " << nu_ << std::endl;
}

// -----------------------------
// Simulation control
// -----------------------------

void MujocoSim::step() {
  mj_step(m_, d_);
}

void MujocoSim::reset() {
  mj_resetData(m_, d_);
  mj_forward(m_, d_);
}

// -----------------------------
// Controller interface
// -----------------------------

ControllerInput MujocoSim::getControllerInput() const {
  ControllerInput in;

  // Floating base:
  // qpos = [7 base | joints]
  // qvel = [6 base | joints]

  const int nq_j = nq_ - 7;
  const int nv_j = nv_ - 6;

  in.q.resize(nq_j);
  in.dq.resize(nv_j);

  in.q    = Eigen::Map<const Eigen::VectorXd>(d_->qpos + 7, nq_j);
  in.dq = Eigen::Map<const Eigen::VectorXd>(d_->qvel + 6, nv_j);

  return in;
}

void MujocoSim::applyTorques(const Eigen::VectorXd& tau) {
  if (tau.size() != nu_) {
    throw std::runtime_error(
      "applyTorques(): tau dimension mismatch"
    );
  }

  Eigen::Map<Eigen::VectorXd>(d_->ctrl, nu_) = tau;
}

// -----------------------------
// Gravity compensation
// -----------------------------

Eigen::VectorXd MujocoSim::computeGravityTorques() {
  // Compute bias forces (gravity + Coriolis)
  mj_inverse(m_, d_);

  // qfrc_bias is size nv (base + joints)
  const int nv_j = nv_ - 6;

  Eigen::VectorXd tau(nv_j);
  tau = Eigen::Map<Eigen::VectorXd>(d_->qfrc_inverse + 6, nv_j);

  return tau;
}
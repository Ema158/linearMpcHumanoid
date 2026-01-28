#include "simulators/MujocoSim.hpp"

#include <iostream>
#include <stdexcept>

// -----------------------------
// Constructor / Destructor
// -----------------------------

MujocoSim::MujocoSim(const std::string& model_path, const Eigen::VectorXd& q0) {
  loadModel(model_path, q0);
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

void MujocoSim::loadModel(const std::string& model_path, const Eigen::VectorXd& q0) {
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

  //Initial position
  auto set_qpos0 = [&](const char* name, double value) {
  int jid = mj_name2id(m_, mjOBJ_JOINT, name);
  if (jid < 0) {
    throw std::runtime_error(std::string("Joint not found: ") + name);
  }
  d_->qpos[m_->jnt_qposadr[jid]] = value;
  };

  // Right leg
  set_qpos0("RHipYawPitch",  q0(6));
  set_qpos0("RHipRoll",      q0(7));
  set_qpos0("RHipPitch",     q0(8));
  set_qpos0("RKneePitch",    q0(9));
  set_qpos0("RAnklePitch",   q0(10));
  set_qpos0("RAnkleRoll",    q0(11));

  // Left leg
  set_qpos0("LHipYawPitch",  q0(12));
  set_qpos0("LHipRoll",      q0(13));
  set_qpos0("LHipPitch",     q0(14));
  set_qpos0("LKneePitch",    q0(15));
  set_qpos0("LAnklePitch",   q0(16));
  set_qpos0("LAnkleRoll",    q0(17));

  //Right arm
  set_qpos0("RShoulderPitch", q0(18));

  //Left arm
  set_qpos0("LShoulderPitch", -q0(23));

  d_->qpos[0] = q0(0);   // x
  d_->qpos[1] = q0(1);   // y
  d_->qpos[2] = q0(2) + 0.035;   // z (above ground) 0.035 cause there is an offset between my base frame and .xml frame

  // base orientation (quaternion)
  d_->qpos[3] = 1.0;   // w
  d_->qpos[4] = 0.0;   // x
  d_->qpos[5] = 0.0;   // y
  d_->qpos[6] = 0.0;   // z
  
  // APPLY IT
  mj_forward(m_, d_);

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
  //std::getchar();
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

void MujocoSim::applyJointPositions(const Eigen::VectorXd& q_des)
{
  kp_.setConstant(nu_, 50.0);
  kd_.setConstant(nu_, 5.0);

  if (q_des.size() != nu_) {
    throw std::runtime_error(
      "applyJointPositions(): q_des dimension mismatch"
    );
  }

  Eigen::VectorXd tau(nu_);
  tau.setZero();

  for (int i = 0; i < nu_; ++i) {
    int joint_id = m_->actuator_trnid[2*i];   // actuator → joint
    int qpos_adr = m_->jnt_qposadr[joint_id];
    int qvel_adr = m_->jnt_dofadr[joint_id];

    double q  = d_->qpos[qpos_adr];
    double dq = d_->qvel[qvel_adr];

    double kp = kp_[i];   // your gains
    double kd = kd_[i];

    tau[i] = kp * (q_des[i] - q);// - kd * dq;
  }

  applyTorques(tau);
}
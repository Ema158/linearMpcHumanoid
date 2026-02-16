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
  set_qpos0("RShoulderRoll", q0(19));
  set_qpos0("RElbowYaw", q0(20));
  set_qpos0("RElbowRoll", q0(21));
  set_qpos0("RWristYaw", q0(22));

  //Left arm
  set_qpos0("LShoulderPitch", -q0(23));
  set_qpos0("LShoulderRoll", q0(24));
  set_qpos0("LElbowYaw", q0(25));
  set_qpos0("LElbowRoll", q0(26));
  set_qpos0("LWristYaw", q0(27));

  //Head
  set_qpos0("HeadYaw", q0(28));
  set_qpos0("HeadPitch", q0(29));

  d_->qpos[0] = q0(0);   // x
  d_->qpos[1] = q0(1);   // y
  d_->qpos[2] = q0(2) + 0.035;   // z (above ground) 0.035 cause there is an offset between my base frame and .xml frame

  // base orientation (quaternion)
  //Curently the orientation of the base is represented as eulerAngles
  Eigen::Vector3d eulerAngles = q0.segment(3,3);
  Eigen::VectorXd quat = eulerAnglesToQuaternion(eulerAngles);
  d_->qpos[3] = quat(0);   // w
  d_->qpos[4] = quat(1);   // x
  d_->qpos[5] = quat(2);   // y
  d_->qpos[6] = quat(3);   // z
  
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

void MujocoSim::applyTorquesV2(const std::vector<jointsIndex>& joints, const Eigen::VectorXd& tau) {
  if (tau.size() != joints.size()) {
    throw std::runtime_error(
      "applyTorques(): tau dimension mismatch"
    );
  }
  
  for(int i = 0; i<joints.size(); i++){
    int index = static_cast<int>(joints[i]);

    if (index < 0 || index >= m_->nu)
            throw std::runtime_error("Invalid actuator index");

    d_->ctrl[index] = tau(i);
  }
}

void MujocoSim::applyJointPositions(const Eigen::VectorXd& q_des)
{
  kp_.setConstant(nu_, 50.0); //50
  kd_.setConstant(nu_, 0.0);

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
    if(i==12){
      std::cout<<"Desired q " << q_des[i] << std::endl;
      std::cout<<"Current q " << q << std::endl;
      std::cout<<"Torque applied = " << tau(i) << std::endl;
    }    
  }
  applyTorques(tau);
}

void MujocoSim::applyJointPositionsV2(const std::vector<jointsIndex>& joints, const Eigen::VectorXd& q_des)
{
  const double kp = 50.0;
  const double kd = 5.0;

  if (q_des.size() != joints.size()) {
    throw std::runtime_error(
      "applyJointPositions(): q_des dimension mismatch"
    );
  }

  Eigen::VectorXd tau(joints.size());
  tau.setZero();

  for (int i = 0; i < joints.size(); ++i) {
    int index = static_cast<int>(joints[i]);

    int joint_id = m_->actuator_trnid[2*index];   // actuator → joint
    int qpos_adr = m_->jnt_qposadr[joint_id];
    int qvel_adr = m_->jnt_dofadr[joint_id];

    double q  = d_->qpos[qpos_adr];
    double dq = d_->qvel[qvel_adr];

    tau[i] = kp * (q_des[i] - q);// - kd * dq;
 
  }
  applyTorquesV2(joints, tau);
}

void MujocoSim::getMujocoState(Eigen::VectorXd& q, Eigen::VectorXd& v)
{
  q = Eigen::Map<Eigen::VectorXd>(d_->qpos, m_->nq);
  v = Eigen::Map<Eigen::VectorXd>(d_->qvel, m_->nv);
}

Eigen::Vector3d MujocoSim::getCoM()
{
  Eigen::Vector3d com;
  com << d_->subtree_com[0],
       d_->subtree_com[1],
       d_->subtree_com[2];
  return com;
}


    
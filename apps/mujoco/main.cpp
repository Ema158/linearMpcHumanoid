#include "simulators/MujocoSim.hpp"
#include "simulators/MujocoViewer.hpp"
#include <iostream>

int main() {
  MujocoSim sim("../models/nao.xml");
  MujocoViewer viewer(sim);

  viewer.run([&]() {
        Eigen::VectorXd tau = sim.computeGravityTorques();
        sim.applyTorques(tau);
    });

  std::cout << "Simulation finished\n";
  return 0;
}
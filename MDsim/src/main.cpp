/// g++ main.cpp MDsim.cpp utils.cpp -o MDsim -O3
/// ./MDsim

#include <cstdlib>   // srand, rand
#include <ctime>     // clock
#include <fstream>   // file
#include <iomanip>   // setprecision
#include <iostream>  // cout

#include "MDSim.hpp"

using namespace std;

void run_sim(string xyz_file = "xyz.in", string run_file = "run.in",
             string thermo_output_file = "thermo.out",
             string traj_output_file = "traj.out") {
  const int Ns = 100;  // output frequency

  MDSim sim_system;
  sim_system.read_run(run_file);
  sim_system.print_run_info();
  sim_system.read_xyz(xyz_file);
  sim_system.print_state();
  sim_system.initializeVelocity(sim_system.temperature);

  const clock_t tStart = clock();
  ofstream thermo_ofile(thermo_output_file);
  thermo_ofile << fixed << setprecision(16);

  ofstream traj_ofile(traj_output_file);
  traj_ofile << fixed << setprecision(16);

  for (size_t step = 0; step < sim_system.numSteps; ++step) {
    if (sim_system.neighbor_flag != 0) {
      sim_system.findNeighbor();
    }

    sim_system.integrate(true, sim_system.timeStep);
    sim_system.findForce();
    sim_system.integrate(false, sim_system.timeStep);
    if (step % Ns == 0) {
      double kineticEnergy = sim_system.findKineticEnergy();
      const double T = kineticEnergy / (1.5 * K_B * sim_system.numAtoms);
      cout << "Step " << step << " T = " << T << " KE = " << kineticEnergy
           << " PE = " << sim_system.pe << endl;
      thermo_ofile << T << " " << kineticEnergy << " " << sim_system.pe << endl;

      sim_system.dump_one_step(traj_ofile, step);
    }
  }

  thermo_ofile.close();
  const clock_t tStop = clock();
  const float tElapsed = float(tStop - tStart) / CLOCKS_PER_SEC;
  cout << sim_system.numUpdates << " neighbor list updates" << endl;
  cout << "Time used = " << tElapsed << " s" << endl;
}

int main(int argc, char *argv[]) {
  /// ./MDsim xyz.in run.in output.out
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " xyz.in run.in" << endl;
    return 1;
  }
  string xyz_file = argv[1];
  string run_file = argv[2];
  string thermo_output_file = argv[3];
  cout << "Running simulation with " << xyz_file << " and " << run_file << endl;
  cout << "Output will be written to " << thermo_output_file << endl;

  run_sim(xyz_file, run_file, thermo_output_file);
  return 0;
}
/// g++ main.cpp MDsim.cpp utils.cpp -o MDsim -O3
/// g++ main.cpp MDsim.cpp utils.cpp -I../include -o MDsim -O3
/// ./MDsim xyz.in run.in thermo.out traj.out

#include <argparse.hpp>
#include <cstdlib>   // srand, rand
#include <ctime>     // clock
#include <fstream>   // file
#include <iomanip>   // setprecision
#include <iostream>  // cout

#include "MDSim.hpp"

using namespace std;

void run_sim(string xyz_file = "xyz.in", string run_file = "run.in",
             string thermo_output_file = "thermo.out",
             string traj_output_file = "traj.out",
             bool dump_xyz_lammps = false) {
  const int Ns = 100;  // output frequency

  MDSim sim_system;
  sim_system.read_run(run_file);
  sim_system.print_run_info();
  sim_system.read_xyz(xyz_file);
  sim_system.print_state();

  sim_system.initializeVelocity(sim_system.temperature);

  if (dump_xyz_lammps) {
    ofstream xyz_lammps_ofile("xyz_lammps.data");
    sim_system.dump_xyz_lammps(xyz_lammps_ofile, 0);
  }
  const clock_t tStart = clock();
  ofstream thermo_ofile(thermo_output_file);
  thermo_ofile << fixed << setprecision(16);

  ofstream traj_ofile(traj_output_file);
  traj_ofile << fixed << setprecision(16);


  // thermo_ofile << "Step Temp(K) Pressure(bar) KE(eV) PE(eV)" << endl;
  for (size_t step = 0; step < sim_system.numSteps; ++step) {
    if (sim_system.neighbor_flag != 0) {
      sim_system.findNeighbor();
    }

    sim_system.integrate(true, sim_system.timeStep);
    sim_system.findForce();
    sim_system.integrate(false, sim_system.timeStep);
    if (step % Ns == 0) {
      sim_system.dump_thermo(cout, step);
      sim_system.dump_thermo(thermo_ofile, step);
      sim_system.dump_trj_one_step(traj_ofile, step);
    }
  }

  thermo_ofile.close();
  const clock_t tStop = clock();
  const float tElapsed = float(tStop - tStart) / CLOCKS_PER_SEC;
  cout << sim_system.numUpdates << " neighbor list updates" << endl;
  cout << "Time used = " << tElapsed << " s" << endl;
}

int main(int argc, char* argv[]) {
  argparse::ArgumentParser program("MDsim");

  program.add_argument("xyz_file").help("Path to the XYZ input file");
  program.add_argument("run_file").help("Path to the run input file");
  program.add_argument("thermo_output_file")
      .help("Path to the thermo output file");
  program.add_argument("traj_output_file")
      .help("Path to the trajectory output file");

  try {
    program.parse_args(argc, argv);
  } catch (const exception& err) {
    cerr << err.what() << endl;
    cerr << program;
    return 1;
  }

  string xyz_file = program.get<string>("xyz_file");
  string run_file = program.get<string>("run_file");
  string thermo_output_file = program.get<string>("thermo_output_file");
  string traj_output_file = program.get<string>("traj_output_file");

  cout << "Running simulation with " << xyz_file << " and " << run_file << endl;
  cout << "Output will be written to " << thermo_output_file << endl;
  cout << "Trajectory will be written to " << traj_output_file << endl;

  run_sim(xyz_file, run_file, thermo_output_file, traj_output_file);
  return 0;
}

#include "MDSim.hpp"

#include <cmath>     // sqrt() function
#include <ctime>     // for timing
#include <fstream>   // file
#include <iomanip>   // setprecision
#include <iostream>  // input/output
#include <iterator>
#include <sstream>  // istringstream
#include <string>   // string
#include <vector>   // vector

using namespace std;

MDSim::MDSim(string run_file, string xyz_file) {
  read_run(run_file);
  read_xyz(xyz_file);
}

void MDSim::read_xyz(string xyz_file) {
  ifstream input(xyz_file);
  if (!input.is_open()) {
    cout << "Failed to open xyz.in." << endl;
    exit(1);
  }

  vector<string> tokens = getTokens(input);

  // line 1
  if (tokens.size() != 1) {
    cout << "The first line of xyz.in should have one item." << endl;
    exit(1);
  }
  numAtoms = getSizeT(tokens[0]);
  cout << "Number of atoms = " << numAtoms << endl;

  // allocate memory
  NN.resize(numAtoms, 0);
  NL.resize(numAtoms * MN, 0);
  mass.resize(numAtoms, 0.0);
  x0.resize(numAtoms, 0.0);
  y0.resize(numAtoms, 0.0);
  z0.resize(numAtoms, 0.0);
  x.resize(numAtoms, 0.0);
  y.resize(numAtoms, 0.0);
  z.resize(numAtoms, 0.0);
  vx.resize(numAtoms, 0.0);
  vy.resize(numAtoms, 0.0);
  vz.resize(numAtoms, 0.0);
  fx.resize(numAtoms, 0.0);
  fy.resize(numAtoms, 0.0);
  fz.resize(numAtoms, 0.0);

  // line 2
  tokens = getTokens(input);
  if (tokens.size() != 9) {
    cout << "The second line of xyz.in should have 9 items." << endl;
    exit(1);
  }

  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = 0; d2 < 3; ++d2) {
      box[d2 * 3 + d1] = getDouble(tokens[d1 * 3 + d2]);
    }
  }
  getInverseBox(box);

  cout << "box matrix H = " << endl;
  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = 0; d2 < 3; ++d2) {
      cout << box[d1 * 3 + d2] << " ";
    }
    cout << endl;
  }

  cout << "inverse box matrix G = " << endl;
  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = 0; d2 < 3; ++d2) {
      cout << box[9 + d1 * 3 + d2] << " ";
    }
    cout << endl;
  }

  // starting from line 3
  for (size_t n = 0; n < numAtoms; ++n) {
    tokens = getTokens(input);
    if (tokens.size() != 5) {
      cout << "The 3rd line and later of xyz.in "
              "should have 5 items."
           << endl;
      exit(1);
    }
    // atom types not used
    x[n] = getDouble(tokens[1]);
    y[n] = getDouble(tokens[2]);
    z[n] = getDouble(tokens[3]);
    mass[n] = getDouble(tokens[4]);
  }

  input.close();
}

void MDSim::read_run(string runfile) {
  ifstream input(runfile);
  if (!input.is_open()) {
    cout << "Failed to open run.in." << endl;
    exit(1);
  }

  while (input.peek() != EOF) {
    vector<string> tokens = getTokens(input);
    if (tokens.size() > 0) {
      if (tokens[0] == "time_step") {
        timeStep = getDouble(tokens[1]);
        if (timeStep < 0) {
          cout << "timeStep should >= 0." << endl;
          exit(1);
        }
        cout << "timeStep = " << timeStep << " fs." << endl;
      } else if (tokens[0] == "run") {
        numSteps = getSizeT(tokens[1]);
        if (numSteps < 1) {
          cout << "numSteps should >= 1." << endl;
          exit(1);
        }
        cout << "numSteps = " << numSteps << endl;
      } else if (tokens[0] == "velocity") {
        temperature = getDouble(tokens[1]);
        if (temperature < 0) {
          cout << "temperature >= 0." << endl;
          exit(1);
        }
        cout << "temperature = " << temperature << " K." << endl;
      } else if (tokens[0] == "neighbor_flag") {
        neighbor_flag = getDouble(tokens[1]);
        if ((neighbor_flag < 0) | (neighbor_flag > 2)) {
          cout << "neighbor_flag can only be 0 or 1 or 2." << endl;
          exit(1);
        }
        cout << "neighbor_flag = " << neighbor_flag << endl;
      } else if (tokens[0][0] != '#') {
        cout << tokens[0] << " is not a valid keyword." << endl;
        exit(1);
      }
    }
  }

  timeStep /= TIME_UNIT_CONVERSION;
  cout << "timeStep = " << timeStep << " natural." << endl;

  input.close();
}

void MDSim::print_run_info() {
  cout << "print run info:" << endl;
  cout << "timeStep = " << timeStep << " natural." << endl;
  cout << "  timeStep = " << timeStep * TIME_UNIT_CONVERSION << " fs." << endl;
  cout << "numSteps = " << numSteps << endl;
  cout << "temperature = " << temperature << " K." << endl;
}

/*
class MDSim {
 public:
  int number;
  int numUpdates = 0;
  int neighbor_flag = 2;
  const int MN = 1000;
  double cutoffNeighbor = 10.0;
  double box[18];
  double pe;
  vector<int> NN, NL;
  vector<double> mass, x0, y0, z0, x, y, z, vx, vy, vz, fx, fy, fz;
}
*/
void MDSim::print_state() {
  cout << "print Sim State:" << endl;
  cout << "number = " << numAtoms << endl;
  cout << "numUpdates = " << numUpdates << endl;
  cout << "neighbor_flag = " << neighbor_flag << endl;
  cout << "MN = " << MN << endl;
  cout << "cutoffNeighbor = " << cutoffNeighbor << endl;
  cout << "box = " << endl;
  for (size_t i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << box[i * 3 + j] << " ";
    }
    cout << endl;
  }
  cout << "pe = " << pe << endl;
  cout << "size of NN = " << NN.size() << endl;
  cout << "size of NL = " << NL.size() << endl;
  cout << "size of mass = " << mass.size() << endl;
}

void MDSim::initializeVelocity(const double T0) {
#ifndef DEBUG
  srand(time(NULL));
#endif
  double centerOfMassVelocity[3] = {0.0, 0.0, 0.0};
  double totalMass = 0.0;
  for (size_t n = 0; n < numAtoms; ++n) {
    totalMass += mass[n];
    vx[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    vy[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    vz[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    centerOfMassVelocity[0] += mass[n] * vx[n];
    centerOfMassVelocity[1] += mass[n] * vy[n];
    centerOfMassVelocity[2] += mass[n] * vz[n];
  }
  centerOfMassVelocity[0] /= totalMass;
  centerOfMassVelocity[1] /= totalMass;
  centerOfMassVelocity[2] /= totalMass;
  for (size_t n = 0; n < numAtoms; ++n) {
    vx[n] -= centerOfMassVelocity[0];
    vy[n] -= centerOfMassVelocity[1];
    vz[n] -= centerOfMassVelocity[2];
  }
  scaleVelocity(T0);
}

void MDSim::findNeighbor() {
  if (checkIfNeedUpdate()) {
    numUpdates++;
    applyPbc();
    if (neighbor_flag == 1)
      findNeighborON1();
    else if (neighbor_flag == 2)
      findNeighborON2();
    updateXyz0();
  }
}

bool MDSim::checkIfNeedUpdate() const {
  bool needUpdate = false;
  for (size_t n = 0; n < numAtoms; ++n) {
    double dx = x[n] - x0[n];
    double dy = y[n] - y0[n];
    double dz = z[n] - z0[n];
    if (dx * dx + dy * dy + dz * dz > 0.25) {
      needUpdate = true;
      break;
    }
  }
  return needUpdate;
}
void MDSim::applyPbc() {
  for (size_t n = 0; n < numAtoms; ++n) {
    // 原始坐标 -> 分数坐标 （相对坐标）
    double sx = box[9] * x[n] + box[10] * y[n] + box[11] * z[n];
    double sy = box[12] * x[n] + box[13] * y[n] + box[14] * z[n];
    double sz = box[15] * x[n] + box[16] * y[n] + box[17] * z[n];
    // 实施 PBC
    applyPbcOne(sx);
    applyPbcOne(sy);
    applyPbcOne(sz);
    // 分数坐标 -> 原始坐标
    x[n] = box[0] * sx + box[1] * sy + box[2] * sz;
    y[n] = box[3] * sx + box[4] * sy + box[5] * sz;
    z[n] = box[6] * sx + box[7] * sy + box[8] * sz;
  }
}

double MDSim::findKineticEnergy() const {
  double kineticEnergy = 0.0;
  for (size_t n = 0; n < numAtoms; ++n) {
    double v2 = vx[n] * vx[n] + vy[n] * vy[n] + vz[n] * vz[n];
    kineticEnergy += mass[n] * v2;
  }
  return kineticEnergy * 0.5;
}

void MDSim::scaleVelocity(const double T0) {
  const double temperature = findKineticEnergy() * 2.0 / (3.0 * K_B * numAtoms);
  double scaleFactor = sqrt(T0 / temperature);
  for (size_t n = 0; n < numAtoms; ++n) {
    vx[n] *= scaleFactor;
    vy[n] *= scaleFactor;
    vz[n] *= scaleFactor;
  }
}

void MDSim::findNeighborON2() {
  const double cutoffSquare = cutoffNeighbor * cutoffNeighbor;
  std::fill(NN.begin(), NN.end(), 0);

  for (size_t i = 0; i < numAtoms - 1; ++i) {
    const double x1 = x[i];
    const double y1 = y[i];
    const double z1 = z[i];
    for (size_t j = i + 1; j < numAtoms; ++j) {
      double xij = x[j] - x1;
      double yij = y[j] - y1;
      double zij = z[j] - z1;
      applyMic(box, xij, yij, zij);
      const double distanceSquare = xij * xij + yij * yij + zij * zij;
      if (distanceSquare < cutoffSquare) {
        NL[i * MN + NN[i]++] = j;
        if (NN[i] > MN) {
          std::cout << "Error: number of neighbors for atom " << i
                    << " exceeds " << MN << std::endl;
          exit(1);
        }
      }
    }
  }
}

void MDSim::findNeighborON1() {
  const double cutoffInverse = 1.0 / cutoffNeighbor;
  double cutoffSquare = cutoffNeighbor * cutoffNeighbor;
  double thickness[3];  // 盒子厚度
  getThickness(thickness, box);

  int numCells[4];

  for (int d = 0; d < 3; ++d) {
    numCells[d] = floor(thickness[d] * cutoffInverse);
  }

  numCells[3] = numCells[0] * numCells[1] * numCells[2];
  int cell[4];

  std::vector<int> cellCount(numCells[3], 0);     // 每个小盒子中的原子数
  std::vector<int> cellCountSum(numCells[3], 0);  // 前缀和

  for (size_t n = 0; n < numAtoms; ++n) {
    const double r[3] = {x[n], y[n], z[n]};
    findCell(box, thickness, r, cutoffInverse, numCells, cell);
    ++cellCount[cell[3]];
  }

  for (int i = 1; i < numCells[3]; ++i) {
    cellCountSum[i] = cellCountSum[i - 1] + cellCount[i - 1];
  }

  std::fill(cellCount.begin(), cellCount.end(), 0);

  std::vector<int> cellContents(numAtoms, 0);

  for (size_t n = 0; n < numAtoms; ++n) {
    const double r[3] = {x[n], y[n], z[n]};
    findCell(box, thickness, r, cutoffInverse, numCells, cell);
    cellContents[cellCountSum[cell[3]] + cellCount[cell[3]]] = n;
    ++cellCount[cell[3]];
  }

  std::fill(NN.begin(), NN.end(), 0);
  // 遍历所有原子，计算作用力时，只考虑近邻列表中的原子
  for (size_t n1 = 0; n1 < numAtoms; ++n1) {
    const double r1[3] = {x[n1], y[n1], z[n1]};
    findCell(box, thickness, r1, cutoffInverse, numCells, cell);
    for (int k = -1; k <= 1; ++k) {
      for (int j = -1; j <= 1; ++j) {
        for (int i = -1; i <= 1; ++i) {
          int neighborCell = cell[3] + (k * numCells[1] + j) * numCells[0] + i;
          if (cell[0] + i < 0) neighborCell += numCells[0];
          if (cell[0] + i >= numCells[0]) neighborCell -= numCells[0];
          if (cell[1] + j < 0) neighborCell += numCells[1] * numCells[0];
          if (cell[1] + j >= numCells[1])
            neighborCell -= numCells[1] * numCells[0];
          if (cell[2] + k < 0) neighborCell += numCells[3];
          if (cell[2] + k >= numCells[2]) neighborCell -= numCells[3];

          for (int m = 0; m < cellCount[neighborCell]; ++m) {
            const size_t n2 = cellContents[cellCountSum[neighborCell] + m];
            if (n1 < n2) {
              double x12 = x[n2] - r1[0];
              double y12 = y[n2] - r1[1];
              double z12 = z[n2] - r1[2];
              applyMic(box, x12, y12, z12);
              const double d2 = x12 * x12 + y12 * y12 + z12 * z12;
              if (d2 < cutoffSquare) {
                NL[n1 * MN + NN[n1]++] = n2;
                if (NN[n1] > MN) {
                  std::cout << "Error: number of neighbors "
                               "for atom "
                            << n1 << " exceeds " << MN << std::endl;
                  exit(1);
                }
              }
            }
          }
        }
      }
    }
  }
}

void MDSim::updateXyz0() {
  for (size_t n = 0; n < numAtoms; ++n) {
    x0[n] = x[n];
    y0[n] = y[n];
    z0[n] = z[n];
  }
}

void MDSim::integrate(const bool isStepOne, const double timeStep) {
  const double timeStepHalf = timeStep * 0.5;
  for (size_t n = 0; n < numAtoms; ++n) {
    const double mass_inv = 1.0 / mass[n];
    const double ax = fx[n] * mass_inv;
    const double ay = fy[n] * mass_inv;
    const double az = fz[n] * mass_inv;
    vx[n] += ax * timeStepHalf;
    vy[n] += ay * timeStepHalf;
    vz[n] += az * timeStepHalf;
    if (isStepOne) {
      x[n] += vx[n] * timeStep;
      y[n] += vy[n] * timeStep;
      z[n] += vz[n] * timeStep;
    }
  }
}

void MDSim::findForce() {
  const double epsilon = 1.032e-2;
  const double sigma = 3.405;
  const double cutoff = 9.0;
  const double cutoffSquare = cutoff * cutoff;
  const double sigma3 = sigma * sigma * sigma;
  const double sigma6 = sigma3 * sigma3;
  const double sigma12 = sigma6 * sigma6;
  const double e24s6 = 24.0 * epsilon * sigma6;
  const double e48s12 = 48.0 * epsilon * sigma12;
  const double e4s6 = 4.0 * epsilon * sigma6;
  const double e4s12 = 4.0 * epsilon * sigma12;
  pe = 0.0;
  for (size_t n = 0; n < numAtoms; ++n) {
    fx[n] = fy[n] = fz[n] = 0.0;
  }

  for (size_t i = 0; i < numAtoms; ++i) {
    const double xi = x[i];
    const double yi = y[i];
    const double zi = z[i];

    if (neighbor_flag == 0) {
      for (size_t j = i + 1; j < numAtoms; ++j) {
        double xij = x[j] - xi;
        double yij = y[j] - yi;
        double zij = z[j] - zi;
        applyMic(box, xij, yij, zij);
        const double r2 = xij * xij + yij * yij + zij * zij;
        if (r2 > cutoffSquare) continue;

        const double r2inv = 1.0 / r2;
        const double r4inv = r2inv * r2inv;
        const double r6inv = r2inv * r4inv;
        const double r8inv = r4inv * r4inv;
        const double r12inv = r4inv * r8inv;
        const double r14inv = r6inv * r8inv;
        const double f_ij = e24s6 * r8inv - e48s12 * r14inv;
        pe += e4s12 * r12inv - e4s6 * r6inv;
        fx[i] += f_ij * xij;
        fx[j] -= f_ij * xij;
        fy[i] += f_ij * yij;
        fy[j] -= f_ij * yij;
        fz[i] += f_ij * zij;
        fz[j] -= f_ij * zij;
      }
    } else {
      for (int jj = 0; jj < NN[i]; ++jj) {
        const int j = NL[i * MN + jj];
        double xij = x[j] - xi;
        double yij = y[j] - yi;
        double zij = z[j] - zi;
        applyMic(box, xij, yij, zij);
        const double r2 = xij * xij + yij * yij + zij * zij;
        if (r2 > cutoffSquare) continue;

        const double r2inv = 1.0 / r2;
        const double r4inv = r2inv * r2inv;
        const double r6inv = r2inv * r4inv;
        const double r8inv = r4inv * r4inv;
        const double r12inv = r4inv * r8inv;
        const double r14inv = r6inv * r8inv;
        const double f_ij = e24s6 * r8inv - e48s12 * r14inv;
        pe += e4s12 * r12inv - e4s6 * r6inv;
        fx[i] += f_ij * xij;
        fx[j] -= f_ij * xij;
        fy[i] += f_ij * yij;
        fy[j] -= f_ij * yij;
        fz[i] += f_ij * zij;
        fz[j] -= f_ij * zij;
      }
    }
  }
}

void MDSim::dump_trj_one_step(std::ostream& os, const size_t step) const {
  os << "ITEM: TIMESTEP\n" << step << "\n";
  os << "ITEM: NUMBER OF ATOMS\n" << numAtoms << "\n";
  os << "ITEM: BOX BOUNDS pp pp pp\n";
  os << box[0] << " " << box[1] << "\n";
  os << box[2] << " " << box[3] << "\n";
  os << box[4] << " " << box[5] << "\n";
  os << "ITEM: ATOMS id type x y z vx vy vz\n";
  size_t atom_type = 1;
  for (size_t n = 0; n < numAtoms; ++n) {
    os << n + 1 << " " << atom_type << " ";
    os << x[n] << " " << y[n] << " " << z[n] << " ";
    os << vx[n] << " " << vy[n] << " " << vz[n] << "\n";
  }
}

void MDSim::dump_xyz_lammps(std::ostream& os, const size_t step) const {
  os << "LAMMPS data file via MDsim, time step " << step << "\n\n";

  os << numAtoms << " atoms\n";
  os << "1 atom types\n\n";
  // box
  // os << box[0] << " " << box[1] << " xlo xhi\n";
  // os << box[2] << " " << box[3] << " ylo yhi\n";
  // os << box[4] << " " << box[5] << " zlo zhi\n\n";
  // double xlo = 0.0, xhi = 10.0, ylo = 0.0, yhi = 10.0, zlo = 0.0, zhi = 10.0;
  // os << xlo << " " << xhi << " xlo xhi\n";
  // os << ylo << " " << yhi << " ylo yhi\n";
  // os << zlo << " " << zhi << " zlo zhi\n\n";
  // masses
  os << "Masses\n\n";
  size_t atom_type = 1;
  double atom_mass = mass[0];
  os << atom_type << " " << atom_mass << "\n\n";
  // atoms
  os << "Atoms\n\n";
  for (size_t i = 1; i <= numAtoms; ++i) {
    os << i << " " << atom_type << " ";
    os << x[i - 1] << " " << y[i - 1] << " " << z[i - 1] << "\n";
  }
  os << "\n";
  // velocities
  os << "Velocities\n\n";
  for (size_t i = 1; i <= numAtoms; ++i) {
    os << i << " " << vx[i - 1] << " " << vy[i - 1] << " " << vz[i - 1] << "\n";
  }
  os << "\n";
}
void MDSim::dump_thermo(std::ostream& os, const size_t step) const {
  double kineticEnergy = findKineticEnergy();
  const double T = kineticEnergy / (1.5 * K_B * numAtoms);
  double pressure = 0.0;
  // T pressure kE pe
  // os << step << " " << T << " " << pressure << " ";
  // os << kineticEnergy << " " << pe << "\n";
  os << T << " " << kineticEnergy << " " << pe << endl;
}

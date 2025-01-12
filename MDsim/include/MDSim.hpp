#pragma once

#include <string>  // string
#include <vector>  // vector
using namespace std;

// Boltzmann's constant in natural unit
const double K_B = 8.617343e-5;
// from natural unit to fs
const double TIME_UNIT_CONVERSION = 1.018051e+1;

class MDSim {
 public:
  /// run members
  int numSteps;
  double temperature;  // K
  double timeStep;     // natural
  /// sim members
  int numAtoms;  // number of atoms
  int numUpdates = 0;
  int neighbor_flag = 2;
  const int maxNeighbors = 1000;  // maximum number of neighbors
  double neighborCutoff = 10.0;   // cutoff of neighbor, in A
  // ax, bx, cx; ay, by, cy; az, bz, cz;
  // remaining 9 elements are the inverse box matrix
  double box[18];
  double potentialEnergy;   // potential energy
  vector<int> numNeighbor;  // neighbor number, neighbor list
  vector<vector<int>>
      neighborLists;  // neighbor list of atom i: neighborLists[i]
  vector<double> mass;
  vector<double> x0, y0, z0;  // initial position, for check delta(r)
  vector<double> x, y, z, vx, vy, vz;
  vector<double> fx, fy, fz;

 private:
  /// @brief 近邻列表线性标度算法
  void findNeighborON1();
  /// @brief 近邻列表平方标度算法
  void findNeighborON2();
  // 周期性边界条件
  void applyPbc();

 public:
  MDSim() = default;
  MDSim(string run_file, string xyz_file);

  void read_xyz(string filename);
  void read_run(string runfile);

  void print_run_info();
  void print_state();

  void initializeVelocity(const double T0);
  void scaleVelocity(const double T0);

  void findNeighbor();
  bool checkIfNeedUpdate() const;

  void updateXyz0();

  void integrate(const bool isStepOne, const double timeStep);
  void findForce();
  double findKineticEnergy() const;

  void dump_thermo(ostream& os, const int step) const;
  void dump_trj_one_step(ostream& os, const int step) const;
  void dump_xyz_lammps(ostream& os, const int step = 0) const;
};

// helper functions
// parse a line of input file and return a vector of tokens
vector<string> getTokens(ifstream& input);
int getInt(string& token);
int getSizeT(string& token);
double getDouble(string& token);

// 计算矩阵 H 的行列式
double getDet(const double box[] /* ax, bx, cx; ay, by, cy; az, bz, cz; */);

// 计算矩阵 H 的逆矩阵
// 逆矩阵的矩阵元对应于数组 box[18] 的后 9 个元素
void getInverseBox(double box[]);

// 计算面积 vector a, b 的面积
float getArea(const double a[3], const double b[3]);
// 求三斜盒子的厚度
void getThickness(double* thickness, double box[18]);

void applyPbcOne(double& sx);

// void findCell(const double box[], const double thickness[], const double r[],
//               double cutoffInverse, const int numCells[], int cell[]);
void findCell(const double box[] /* ax, bx, cx; ay, by, cy; az, bz, cz; ... */,
              const double thickness[] /* h1, h2, h3; */,
              const double r[] /* x, y, z; of the atom */,
              double cutoffInverse /* 1.0 / neighborCutoff*/,
              const int numCells[] /* Na, Nb, Nc, N; */,
              int cell[] /* i, j, k, index; */);

// 相对分数坐标实施最小镜像约定
void applyMicOne(double& x12);

// 三斜盒子的最小镜像约定
void applyMic(const double box[], double& x12, double& y12, double& z12);
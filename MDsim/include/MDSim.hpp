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
  size_t numSteps;
  double temperature;  // K
  double timeStep;     // natural
  /// sim members
  size_t numAtoms;  // number of atoms
  int numUpdates = 0;
  int neighbor_flag = 2;
  const int MN = 1000;  // maximum number of neighbors
  double cutoffNeighbor = 10.0;
  // ax, bx, cx; ay, by, cy; az, bz, cz;
  // remaining 9 elements are the inverse box matrix
  double box[18];
  double pe; // potential energy
  std::vector<int> NN, NL;
  std::vector<double> mass, x0, y0, z0, x, y, z, vx, vy, vz, fx, fy, fz;

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

  void dump_thermo(std::ostream& os, const size_t step) const;
  void dump_trj_one_step(std::ostream& os, const size_t step) const;
  void dump_xyz_lammps(std::ostream& os, const size_t step = 0) const;
};

// helper functions
// parse a line of input file and return a vector of tokens
std::vector<std::string> getTokens(std::ifstream& input);
int getInt(std::string& token);
size_t getSizeT(std::string& token);
double getDouble(std::string& token);

// 计算矩阵 H 的行列式
double getDet(const double box[]);

// 计算矩阵 H 的逆矩阵
// 逆矩阵的矩阵元对应于数组 box[18] 的后 9 个元素
void getInverseBox(double box[]);

// 计算面积 vector a, b 的面积
float getArea(const double a[3], const double b[3]);
// 求三斜盒子的厚度
void getThickness(double* thickness, double box[18]);

void applyPbcOne(double& sx);

void findCell(const double box[], const double thickness[], const double r[],
              double cutoffInverse, const int numCells[], int cell[]);

// 相对分数坐标实施最小镜像约定
void applyMicOne(double& x12);
// 三斜盒子的最小镜像约定
void applyMic(const double box[], double& x12, double& y12, double& z12);
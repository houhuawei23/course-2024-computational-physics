#include <cmath>     // sqrt() function
#include <ctime>     // for timing
#include <fstream>   // file
#include <iomanip>   // std::setprecision
#include <iostream>  // input/output
#include <iterator>
#include <sstream>  // std::istringstream
#include <string>   // string
#include <vector>   // vector
using namespace std;

std::vector<std::string> getTokens(std::ifstream& input) {
  std::string line;
  std::getline(input, line);
  std::istringstream iss(line);
  std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                  std::istream_iterator<std::string>{}};
  return tokens;
}

size_t getSizeT(std::string& token) {
  size_t value = 0;
  try {
    value = std::stoll(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:" << e.what() << std::endl;
    exit(1);
  }
  return value;
}

int getInt(std::string& token) {
  int value = 0;
  try {
    value = std::stoi(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:" << e.what() << std::endl;
    exit(1);
  }
  return value;
}

double getDouble(std::string& token) {
  float value = 0;
  try {
    value = std::stod(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:" << e.what() << std::endl;
    exit(1);
  }
  return value;
}

double getDet(const double box[]) {
  return box[0] * (box[4] * box[8] - box[5] * box[7]) +
         box[1] * (box[5] * box[6] - box[3] * box[8]) +
         box[2] * (box[3] * box[7] - box[4] * box[6]);
}

void getInverseBox(double box[]) {
  box[9] = box[4] * box[8] - box[5] * box[7];
  box[10] = box[2] * box[7] - box[1] * box[8];
  box[11] = box[1] * box[5] - box[2] * box[4];
  box[12] = box[5] * box[6] - box[3] * box[8];
  box[13] = box[0] * box[8] - box[2] * box[6];
  box[14] = box[2] * box[3] - box[0] * box[5];
  box[15] = box[3] * box[7] - box[4] * box[6];
  box[16] = box[1] * box[6] - box[0] * box[7];
  box[17] = box[0] * box[4] - box[1] * box[3];
  double det = getDet(box);
  for (int n = 9; n < 18; ++n) {
    box[n] /= det;
  }
}

float getArea(const double a[3], const double b[3]) {
  const double s1 = a[1] * b[2] - a[2] * b[1];
  const double s2 = a[2] * b[0] - a[0] * b[2];
  const double s3 = a[0] * b[1] - a[1] * b[0];
  return sqrt(s1 * s1 + s2 * s2 + s3 * s3);
}

void getThickness(double* thickness, double box[18]) {
  double volume = abs(getDet(box));
  const double a[3] = {box[0], box[3], box[6]};
  const double b[3] = {box[1], box[4], box[7]};
  const double c[3] = {box[2], box[5], box[8]};
  thickness[0] = volume / getArea(b, c);
  thickness[1] = volume / getArea(c, a);
  thickness[2] = volume / getArea(a, b);
}
/*
确定粒子所在的格子编号
*/
void findCell(const double box[], const double thickness[], const double r[],
              double cutoffInverse, const int numCells[], int cell[]) {
  double s[3];
  s[0] = box[9] * r[0] + box[10] * r[1] + box[11] * r[2];
  s[1] = box[12] * r[0] + box[13] * r[1] + box[14] * r[2];
  s[2] = box[15] * r[0] + box[16] * r[1] + box[17] * r[2];
  for (int d = 0; d < 3; ++d) {
    cell[d] = floor(s[d] * thickness[d] * cutoffInverse);
    if (cell[d] < 0) cell[d] += numCells[d];
    if (cell[d] >= numCells[d]) cell[d] -= numCells[d];
  }
  cell[3] = cell[0] + numCells[0] * (cell[1] + numCells[1] * cell[2]);
}

void applyPbcOne(double& sx) {
  if (sx < 0.0) {
    sx += 1.0;
  } else if (sx > 1.0) {
    sx -= 1.0;
  }
}

void applyMicOne(double& x12) {
  if (x12 < -0.5)
    x12 += 1.0;
  else if (x12 > +0.5)
    x12 -= 1.0;
}

void applyMic(const double box[], double& x12, double& y12, double& z12) {
  // 原始坐标 -> 分数坐标 （相对坐标）
  double sx12 = box[9] * x12 + box[10] * y12 + box[11] * z12;
  double sy12 = box[12] * x12 + box[13] * y12 + box[14] * z12;
  double sz12 = box[15] * x12 + box[16] * y12 + box[17] * z12;
  // 标实施最小镜像约定
  applyMicOne(sx12);
  applyMicOne(sy12);
  applyMicOne(sz12);
  // 分数坐标 -> 原始坐标
  x12 = box[0] * sx12 + box[1] * sy12 + box[2] * sz12;
  y12 = box[3] * sx12 + box[4] * sy12 + box[5] * sz12;
  z12 = box[6] * sx12 + box[7] * sy12 + box[8] * sz12;
}

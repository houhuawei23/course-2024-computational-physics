# [course 2024] Computational Physics 

reference:

- [Molecular-Dynamics-Simulation: 樊哲勇](https://github.com/brucefan1983/Molecular-Dynamics-Simulation)
- [github.com: p-ranav/argparse](https://github.com/p-ranav/argparse)

## MDsim

> input & init -> sim n steps -> post-processing

build

```bash
cd MDsim
make
./bin/MDsim xyz.in run.in thermo.out traj.out

cd ./test/ArgonCrystal/
../../bin/MDsim xyz.in run.in therom.out traj.out
```

## LammpsScripts

Scripts for Lammps simulation.

- `LammpsModeling/` Lammps 建模练习
- `BrownMotion/` 布朗运动
  - Class09, 布朗运动 计算扩散系数
  - 热浴 - 用于控温
  - Berendsen 热浴、Bussi-Donadio-Parrinello 热浴、
  - Nose-Hoover 热浴、Nose-Hoover 链热浴以及朗之万热浴
- `Diffusion/` 扩散运动
  - Class09, 布朗运动 计算扩散系数
  - Ar的扩散系数
- `VacuumDiffusion/` 真空扩散
- `Flow/` 流动 Class10
  - 流速随距离 x 的变化
- `Vicosity/` 粘度 Class10
- `BalanceLatticeConstant/` 平衡晶格常数
  - Class11, 力学性质计算
  - 循环控制 二次曲线拟合
- `LJ-BodyModulus/` LJ 体模量
  - Class11, 力学性质计算
  - BM (Birch-Murnaghan) 方程拟合
- `ThermalExpansion/` 热膨胀
  - Class12, 热力学性质计算
- `MeltingPoint/` 熔点
  - Class12
- `HeatConduction/` 热传导
  - Class13, 热输运性质计算
  - 傅里叶定律; 非平衡法计算热导率; Green-Kubo理论 


others

- `potentials/` 相互作用势

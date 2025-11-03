# Unified Library for Kinematic Identification and Inverse Kinematics in C++ and Matlab

**DMH** is a fast and robust kinematics library for serial manipulators. 
It aims to establish a unified framework for kinematic identification (KI) and inverse kinematics (IK) and to develop an efficient, robust, and accurate solution method for both problems.
It is a **fifth-order POE-based** method, **published in IEEE-TRO (conditionally accepted)**, that enhances the identification accuracy for solving KI problems and improves the convergence rate and robustness for solving IK problems.

> Yuhan Chen, et al., "A Fifth-Order POE-Based Method for Kinematic Identification and Inverse Kinematics of Serial Robots," *IEEE Transactions on Robotics* (conditionally accepted).

## Motivation

* Illustrated the damped modified Halley (DMH) method (fifth-order convergence) for KI and IK
* Disseminate published academic results

The MATLAB codebase for the DMH method provides a unified function, named *DMH_Method_KI_IK*, that consistently solves both KI and IK problems, following the procedure described in *Algorithm 1* of the paper.
The C++ codebase implements the DMH, DTH, and DNR methods for solving IK problems.
Qt is recommended for use with the C++ examples, as the project configuration files are provided in the “.pro” format.
**Note**: The InverseKinematics class operates independently of Qt and can be used in standalone C++ applications.

## Prerequisites

* Qt (*[Markdown Guide](https://www.markdownguide.org)*)

## How it works

本项目基于已发表于 TRO 的算法实现，以 DH 建模为基础，利用数值迭代高效求解 IK。实现细节与理论推导请参考论文与 `docs/` 中的材料。

## 如何构建

本项目使用 CMake 进行构建。为了获得最佳性能，请使用 `Release` 模式进行构建。

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## 代码说明与用法

- 头文件与核心接口：
	- `your_project/Robot.hpp`：机器人结构、正向/雅可比计算；
	- `your_project/IKSolver.hpp`：逆运动学求解接口。
- 示例：`src/sample_usage.cpp` 展示最小化用法（创建 Robot 与 IKSolver，给定目标位姿与初值，得到解）。

### 基本用法

本库使用 **Denavit-Hartenberg (DH)** 参数来定义机器人结构。通过提供 DH 参数、关节类型等来定义 `Robot` 和 `IKSolver` 对象。

```cpp
// 1. 包含头文件
#include "your_project/Robot.hpp"
#include "your_project/IKSolver.hpp"

// 2. 定义机器人参数 (DH, 关节类型等)
auto robot = std::make_shared<your_project::Robot<DOF>>(...);

// 3. 创建逆运动学求解器
const your_project::IKSolver<DOF> ik_solver(robot);

// 4. 准备目标位姿和初始猜测角
// ...

// 5. 求解
ik_solver.solve(target_pose, q_initial_guess, q_solution);
```

更详细的用法请参考 `src/sample_usage.cpp` 文件。

## 依赖项

- 线性代数库：**Eigen 3.4**（必需）

如果你的系统尚未安装（Ubuntu）：

```bash
sudo apt install libeigen3-dev
```

## 如何引用

如果你在研究或产品中使用了本项目，请引用我们的论文（发表于 IEEE Transactions on Robotics, TRO）：

[1] [作者], [论文标题], IEEE Transactions on Robotics (TRO), [年份]. doi: [DOI]

## 许可

本项目在 **AGPL-3.0** 许可下发布。若你的商业场景无法开源衍生作品，可邮件联系获取商业授权：`your.email@domain.com`。

## 商业许可

本项目默认遵循 AGPL-3.0。若你的应用不便开源衍生作品（例如闭源分发或仅以服务形式提供），欢迎联系洽谈商业授权：`your.email@domain.com`。

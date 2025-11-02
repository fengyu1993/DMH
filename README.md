# Unified Library for Kinematic Identification and Inverse Kinematics in C++ and Matlab

**DMH** is a fast and robust kinematics library for serial manipulators. 
It aims to establish a unified framework for kinematic identification (KI) and inverse kinematics (IK) and to develop an efficient, robust, and accurate solution method for both problems.
It is based on the novel damped modified Halley method, published in **IEEE Transactions on Robotics (TRO)**, that enhances the identification accuracy for solving KI problems and improves the convergence rate
and robustness for solving IK problems.

> 在此处填写你的论文引用信息。
> 例如：[作者], "[论文标题]," IEEE Transactions on Robotics (TRO), [年份].

## 算法简介

本库实现了基于 **[你的算法名称]** 的广义逆运动学求解器，面向串联机械臂的高效与高鲁棒应用场景。该方法在迭代求解策略与误差建模上进行了改进，相较传统一阶方法可更快收敛、失败率更低，适合实时与高可靠性要求的任务。

## 性能与基准

以下为在典型 6-DOF 机械臂、较大初始关节误差条件下的示例基准（请替换为你的数据）：

| 求解器       | 平均求解时间 | 错误率 |
| ------------ | ------------ | ------ |
| 你的算法     | XX μs        | X.XX%  |
| 对比算法 A   | YY μs        | Y.YY%  |
| 对比算法 B   | ZZ μs        | Z.ZZ%  |

当初值较好（接近真解）时，错误率可趋近于零，求解时间也将显著降低。

更多细节可参考论文与预印本（可选）：`docs/paper_preprint.pdf`。

> 推荐引用：
> [作者], "[论文标题]," IEEE Transactions on Robotics (TRO), [年份]. doi: [DOI]

## 本项目功能

### 本项目提供

- 高效、鲁棒的广义 **逆运动学**（IK）求解（6 自由度工具位姿约束：位置 + 朝向）。
- **正向运动学**（FK）计算，返回末端或任意关节的位姿。
- **速度运动学/雅可比矩阵**（Jacobian）计算，用于速度映射与敏感度分析。
- 支持固定 DOF 模板与动态 DOF 两种构建方式，便于在实时/嵌入式场景下避免动态分配。

### 不包含/暂未提供

- 零空间优化（针对冗余机械臂的二次目标优化）。
- 降维/部分约束 IK（例如仅位置或仅姿态等子任务）。
- 可行性/安全性检查（关节/速度/加速度/工作空间约束）。
- 碰撞检测、轨迹规划与控制功能（超出运动学范畴）。

### 适用范围与前提

- 串联机械臂（单链，不含分支/闭环）。
- 采用 Denavit–Hartenberg（DH）建模（Spong 记号与编号约定）。
- 支持旋转/移动关节混合，支持自定义基座与工具坐标系。
- 建议统一单位（米、弧度）以避免数值缩放问题。

## 为何使用通用逆运动学

如果你的机器人具备解析解（封闭解），优先使用解析解；但在以下情况，通用 IK 更为合适：

- 无法或不便推导解析解；
- 机器人经过标定，DH 参数扰动导致解析解失效；
- 结构复杂、无球腕等，位置与姿态难以解耦；
- 需要后续扩展零空间优化等能力（本仓库计划支持）。

## 工作原理

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

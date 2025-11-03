# Unified Library for Kinematic Identification and Inverse Kinematics in C++ and Matlab

**DMH** is a fast and robust kinematics library for serial manipulators. 
It aims to establish a unified framework for kinematic identification (KI) and inverse kinematics (IK) and to develop an efficient, robust, and accurate solution method for both problems.
It is a fifth-order Product-of-Exponential (POE)-based method, **published in IEEE-TRO (conditionally accepted)**, that enhances the identification accuracy for solving KI problems and improves the convergence rate and robustness for solving IK problems.

> Yuhan Chen, et al., "A Fifth-Order POE-Based Method for Kinematic Identification and Inverse Kinematics of Serial Robots," *IEEE Transactions on Robotics* (conditionally accepted).

## Motivation

* Illustrated the damped modified Halley (DMH) method (fifth-order convergence) for KI and IK
* Disseminate published academic results

The MATLAB codebase for the DMH method provides a unified function, named *DMH_Method_KI_IK*, that consistently solves both KI and IK problems, following the procedure described in *Algorithm 1* of the paper.
The C++ codebase implements the DMH, DTH, and DNR methods for solving IK problems.
Qt is recommended for use with the C++ examples, as the project configuration files are provided in the “.pro” format.
**Note**: The InverseKinematics class operates independently of Qt and can be used in standalone C++ applications.

## Innovation 

The DMH repository is based on the DMH algorithm, the full details of which can be found in the paper docs/SLloydEtAl2022_QuIK_preprint.pdf.

本项目基于已发表于 TRO 的算法实现，以 DH 建模为基础，利用数值迭代高效求解 IK。实现细节与理论推导请参考论文与 `docs/` 中的材料。

## Prerequisites

* Qt (*[https://www.qt.io/download](https://www.qt.io/download)*)
* MATLAB (*[https://www.mathworks.com/products/matlab.html](https://www.mathworks.com/products/matlab.html)*)
* Eigen3 (*[ http://eigen.tuxfamily.org]( http://eigen.tuxfamily.org)*)

A copy of Eigen (version 3.4) is included in the `dependencies/` directory for convenience.

## Contributing

Feel free to submit pull requests and use the issue tracker to start a discussion about any bugs you encounter. Please provide a description of your compiler and operating system for any software related bugs.

## Citations

If you use our work, please reference our publication below. Recommended citation:
> [1] XXXXXXXX. doi: XXXXXXX

## License

This project is licensed under the MIT License - see the *[LICENSE](https://github.com/fengyu1993/DMH/blob/main/LICENSE)* file for details.


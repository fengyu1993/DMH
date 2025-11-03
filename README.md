# Unified Library for Kinematic Identification and Inverse Kinematics in C++ and Matlab

**DMH** is a fast and robust kinematics library for serial manipulators. 
It aims to establish a unified framework for kinematic identification (KI) and inverse kinematics (IK) and to develop an efficient, robust, and accurate solution method for both problems.
It is a fifth-order Product-of-Exponential (POE)-based method, **published in IEEE-TRO (conditionally accepted)**, that enhances the identification accuracy for solving KI problems and improves the convergence rate and robustness for solving IK problems.

> Yuhan Chen, et al., "A Fifth-Order POE-Based Method for Kinematic Identification and Inverse Kinematics of Serial Robots," *IEEE Transactions on Robotics* (conditionally accepted).

## Motivation

* Illustrated the DMH method (fifth-order convergence) for KI and IK
* Disseminate published academic results

The C++ codebase implements the damped modified Halley (DMH), damped traditional Halley (DTH), and damped Newton–Raphson (DNR) methods for solving IK problems.
Qt is recommended for use with the C++ examples, as the project configuration files are provided in the “.pro” format.
The MATLAB codebase for the DMH method provides a unified function, named *DMH_Method_KI_IK*, that consistently solves both KI and IK problems, following the procedure described in *Algorithm 1* of the paper.

**Note**: The InverseKinematics class operates independently of Qt and can be used in standalone C++ applications.

## Innovation 

The DMH repository is based on the DMH algorithm, the full details of which can be found in the paper *[docs/YuhanChen2025_DMH_preprint.pdf](https://github.com/fengyu1993/DMH/blob/main/LICENSE).

The key innovation of this work lies in incorporating higher-order derivatives into the iterative solver, where the use of the Hessian matrix extends the convergence rate from third order to fifth order.
Most existing KI and IK solvers employ the geometric Jacobian of the kinematic chain within a damped Newton framework, where the robot's kinematic function is linearized at each iteration, and this local linear approximation is used to project the estimate toward the desired solution.
The *[DQuIK](https://github.com/steffanlloyd/quik)* method introduces a third-order iterative scheme that leverages not only the Jacobian but also the Hessian matrix.


## Prerequisites

* Qt (*[https://www.qt.io/download](https://www.qt.io/download)*)
* MATLAB (*[https://www.mathworks.com/products/matlab.html](https://www.mathworks.com/products/matlab.html)*)
* Eigen3 (*[ http://eigen.tuxfamily.org]( http://eigen.tuxfamily.org)*)

A copy of Eigen (version 3.4) is included in the `dependencies/` directory for convenience.

## Contributing

Feel free to submit pull requests and use the issue tracker to start a discussion about any bugs you encounter. Please describe your compiler and operating system for any software-related bugs.

## Citations

If you use our work, please reference our publication below. Recommended citation:
> [1] XXXXXXXX. doi: XXXXXXX

## License

This project is licensed under the MIT License - see the *[LICENSE](https://github.com/fengyu1993/DMH/blob/main/LICENSE)* file for details.


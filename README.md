# Unified Library for Kinematic Identification and Inverse Kinematics in C++ and Matlab

**DMH** is a fast and robust kinematics library for serial manipulators. 
It aims to establish a unified framework for kinematic identification (KI) and inverse kinematics (IK) and to develop an efficient, robust, and accurate solution method for both problems.
It is a fifth-order Product-of-Exponential (POE)-based method, **published in IEEE-TRO (conditionally accepted)**, that enhances the identification accuracy for solving KI problems and improves the convergence rate and robustness for solving IK problems.

> Yuhan Chen, et al., "A Fifth-Order POE-Based Method for Kinematic Identification and Inverse Kinematics of Serial Robots," *IEEE Transactions on Robotics* (conditionally accepted).

## Motivation

* Illustrated the DMH method (fifth-order convergence) for KI and IK
* Disseminate published academic results

The C++ codebase compiled with the **MSVC 2019 compiler** implements the damped modified Halley (DMH), damped traditional Halley (DTH), and damped Newton–Raphson (DNR) methods for solving IK problems.
**Qt** is recommended for use with the C++ examples, as the project configuration files are provided in the “.pro” format.
The MATLAB codebase for the DMH method provides a unified function, named *DMH_Method_KI_IK*, that consistently solves both KI and IK problems, following the procedure described in *Algorithm 1* of the paper.

**Note**: The InverseKinematics class operates independently of Qt and can be used in standalone C++ applications.

## Innovation 

The DMH repository is based on the DMH algorithm, the full details of which can be found in the paper *[docs/YuhanChen2025_DMH_preprint.pdf](https://github.com/fengyu1993/DMH/blob/main/LICENSE)*.

The key innovation of this work lies in incorporating higher-order derivatives into the iterative solver, where the use of the Hessian matrix extends the convergence rate from third order to fifth order.
Most existing KI and IK solvers employ the geometric Jacobian of the kinematic chain within a damped Newton framework, where the robot's kinematic function is linearized at each iteration, and this local linear approximation is used to project the estimate toward the desired solution.
The *[DQuIK](https://github.com/steffanlloyd/quik)* (DTH) method introduced a third-order iterative scheme that leverages not only the Jacobian but also the Hessian matrix.
The key point of the DQuIK method lies in its use of first- and second-order derivatives to approximate the robot’s kinematic function with a quadratic model, which is then used to project toward the solution.
Our innovation is to utilize both the Jacobian and Hessian matrices to predict the next-substep derivative based on the DTH update, after which a single DNR step is performed to reach the next result.
It can be observed that the proposed DMH method only adds one evaluation of the error at another point iterated by the DTH methods. Still, its order of convergence will be proven to increase from three to five in the paper.
This is visualized in the image below.

<img width="4002" height="670" alt="图片1" src="https://github.com/user-attachments/assets/fc7662c3-02d4-40cf-85ed-362a5665cd22" />

<div align="center">
(a) 
</div>

<img width="4002" height="662" alt="图片2" src="https://github.com/user-attachments/assets/dff87c80-c8d2-4794-be1e-191d9d8a9301" />

<div align="center">
(b) 
</div>

Fig. Graphical illustration of the NR, TH, and MH methods on one-dimensional functions. 
(a) $\Xi = \arctan(\Phi)$. (b) $\Xi = \sin(\Phi) + 0.5\cos(3\Phi) - 0.1\Phi$.
The MH method achieves the fewest iterations (top), while the NR and TH methods diverge, with only the MH method converging (bottom). In the latter case, the $\Xi'(\Phi_{\mathrm{TH}})$ estimated by $\Xi''(\Phi_0)$ has the opposite sign of the true $\Xi'(\Phi_{\mathrm{TH}})$, which facilitates convergence.

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


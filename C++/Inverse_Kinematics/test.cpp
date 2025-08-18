#include <iostream>
#include "InverseKinematics.h"


void   getSCARAParameter(MatrixXd &Slist, VectorXd &Gamma);
void   getM710iCParameter(MatrixXd &Slist, VectorXd &Gamma);

int main()
{
    // Initialization
    int flag = 1;
    MatrixXd Slist;
    VectorXd Gamma;
    double epsilon = 1e-20;
    int tau_max = 30;
    double sigma = 1e-6;
    srand(1);
    //
    VectorXd delta_q_init, delta_q_desired;
    switch(flag)
    {
    case 1:
        getSCARAParameter(Slist, Gamma);
          break;
    case 2:
         getM710iCParameter(Slist, Gamma);
        break;
    }
    delta_q_init.resize(Slist.cols());
    delta_q_desired.resize(Slist.cols());
    delta_q_init <<  VectorXd::Random(Slist.cols()) * PI;
    delta_q_desired << VectorXd::Random(Slist.cols()) * PI;

    InverseKinematics IK = InverseKinematics(Slist, Gamma, epsilon, tau_max, sigma);
    Matrix4d Td = IK.getTeFK(Slist, IK.Exp6(Gamma), delta_q_desired);


    // DMH
    IK.setAlgorithm(ALG_DMH);
    IK.solver(delta_q_init, Td);
    cout << "algo: " <<  IK.getAlgorithm() << endl;
    cout << "delta_q: \n" <<  IK.getResult() << endl;
    cout << "error: \n" <<  IK.getError() << endl;

    return 0;
}


void   getM710iCParameter(MatrixXd &Slist, VectorXd &Gamma)
{
    Slist.resize(6,6);
    Slist << 0,   0,      0,   1.0,    0,   1.0,
        0,   1.0,   1.0,   0,   1.0,   0,
        1.0,   0,   0,   0,   0,   0,
        0,  -0.565,  -1.435,  -0,  -1.605,  -0,
        -0.5,  -0,  -0,   1.605,  -0,   1.605,
        0,   0.650,   0.650,   0,   1.666045156698660,   0;
    Gamma.resize(6);
    Gamma << -1.209199576156145,   1.209199576156145,  -1.209199576156145,   0.688951829140432,  -0.252376178687183,   2.629717148871045;
}


void   getSCARAParameter(MatrixXd &Slist, VectorXd &Gamma)
{
    Slist.resize(6,4);
    Slist << 0,                   0,                   0,                   0,
                 0,                   0,                   0,                   0,
                 1,                   1,                   0,                    -1,
                 0,                   0,                   0,                   0,
                 0,                 -0.25,               0,                 0.47,
                 0,                   0,                    -1,                   0;
    Gamma.resize(6);
    Gamma << 0,           0,           0,           0,           0,        0;
}

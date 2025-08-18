#ifndef INVERSEKINEMATICS_H
#define INVERSEKINEMATICS_H

#define  EPS_NEAR_ZERO      1e-7
#define  PI          3.14159265358979323846

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

enum AlgorithmType {
    ALG_DMH = 0,
    ALG_DTH,
    ALG_DNR
};

class InverseKinematics
{
public:
    InverseKinematics(const MatrixXd &Slist, const VectorXd &Gamma, double epsilon = 1e-13, int tau_max = 30, double sigma = 1e-6);
    void solver(const VectorXd &delta_q_init, const Matrix4d &Td);

    bool nearZero(double v);
    void TransToRp(const Matrix4d &T, Matrix3d &R, Vector3d &p);
    Matrix4d TransInv(const Matrix4d &T);
    Matrix3d vecToSo3(const Vector3d &v);
    Vector3d so3ToVec(const Matrix3d &m);
    Matrix<double,6,6> adjointInv(const Matrix4d &T);
    Matrix<double,6,6> ad(const Matrix<double,6,1> &V);
    Matrix3d getDexpW(const Vector3d &lambda_w);
    Matrix3d getDexpWInv(const Vector3d &lambda_w);
    Matrix<double,6,6> getDexpInv(const Matrix<double,6,1> &lambda);
    Matrix4d Exp6(const Matrix<double,6,1> &lambda);
    Matrix<double,6,1> Log6(const Matrix4d &T);

    void setAlgorithm(AlgorithmType algo);
    void setDeltaQInit(const VectorXd &delta_q_init);
    void setTd(Matrix4d &Td);
    void setEpsilon(double epsilon);
    void setTauMax(int tau_max);
    void setSigma(double sigma);
    string getAlgorithm();
    VectorXd getDeltaQInit();
    Matrix4d getTd();
    double getEpsilon();
    int getTauMax();
    double getSigma();
    VectorXd getResult();
    VectorXd getError();

    void getTlist(const MatrixXd &Slist, const MatrixXd& M, const VectorXd &delta_q, vector<Matrix4d> &T_list);
    void getTinlist(const vector<Matrix4d> &Tlist,  int n, vector<Matrix4d> &Tinlist);
    void getJacobian(const MatrixXd &Slist, const vector<Matrix4d> &Tin_list, MatrixXd &J);
    void getHessian(const MatrixXd &J, vector<MatrixXd>& Hessian);
    void getLeftInverse(const MatrixXd &A,    const VectorXd &b,    VectorXd& x);
    Matrix4d getTeFK(const MatrixXd &Slist, const MatrixXd& M, const VectorXd &delta_q);

private:
    MatrixXd Slist_;
    VectorXd Gamma_;
    VectorXd delta_q_init_;
    int n_;
    Matrix4d Td_;
    double epsilon_;
    int tau_max_;
    double sigma_;
    vector<Matrix4d> Tlist_;
    MatrixXd I_;
    vector<double> error_;
    AlgorithmType algo_;
    VectorXd delta_q_;
    bool flag_;
    Matrix4d M_;
    vector<Matrix4d> Tinlist_;
    MatrixXd JXi_;
    vector<MatrixXd> Hessian_;
    VectorXd DeltaQDNR_;
    VectorXd DeltaQDTH_;
    VectorXd DeltaQDMH_;
    MatrixXd tempSum_;
};

#endif // INVERSEKINEMATICS_H

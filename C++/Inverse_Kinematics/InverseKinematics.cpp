/*
 * Copyright (c) 2025 Yuhan Chen (chenyuhan19930920@163.com)
 * Licensed under the MIT License.
 * See LICENSE file in the project root for full license information.
*/

#include "InverseKinematics.h"
#include <iostream>

InverseKinematics::InverseKinematics(const MatrixXd &Slist, const VectorXd &Gamma, double epsilon, int tau_max, double sigma)
    :   Slist_(Slist),                  Gamma_(Gamma),             n_((int)Slist.cols()),       epsilon_(epsilon),                   tau_max_(tau_max),
        sigma_(sigma),              Tlist_(this->n_+1, Matrix4d::Identity()),               I_(MatrixXd::Identity(this->n_, this->n_)),        algo_(ALG_DMH),
        delta_q_(VectorXd::Zero(this->n_)),                       flag_(false),                    M_(Exp6(this->Gamma_)),   Tinlist_(this->n_+1, Matrix4d::Identity()),
        JXi_(MatrixXd::Zero(6, this->n_)),                         Hessian_(this->n_, MatrixXd::Zero(6, this->n_)),
        DeltaQDNR_(VectorXd::Zero(this->n_)),               DeltaQDTH_(VectorXd::Zero(this->n_)),                    DeltaQDMH_(VectorXd::Zero(this->n_)),
        tempSum_(MatrixXd::Zero(6, this->n_))
{}

void InverseKinematics::solver(const VectorXd &delta_q_init, const Matrix4d &Td)
{
    this->delta_q_ = delta_q_init;
    this->Td_ = Td;
    this->error_.clear();
    int tau;
    for (tau = 0; tau < this->tau_max_; ++tau)
    {
        // Eq. 16
        getTlist(this->Slist_, this->M_, this->delta_q_, this->Tlist_);
        getTinlist(this->Tlist_, this->n_, this->Tinlist_);
        Matrix4d Te =   this->Tinlist_[0];
        // Eq. 61
        Matrix<double,6,1> xi = Log6(TransInv(Te) * this->Td_);
        double err = xi.dot(xi);
        this->error_.push_back(err);
        if (err < this->epsilon_)
            break;
        // Eq. 87
        getJacobian(this->Slist_, this->Tinlist_, this->JXi_);
        // Methods
        switch (this->algo_)
        {
        case ALG_DMH:
            // Eq. 30
            getLeftInverse(this->JXi_, -xi, this->DeltaQDNR_);
            // Eq. 92
            getHessian(this->JXi_, this->Hessian_);
            // // Eq. 39
            this->tempSum_.setZero();
            for (int i = 0; i < this->n_; ++i) {
                this->tempSum_ += this->DeltaQDNR_(i) * this->Hessian_[i];
            }
            getLeftInverse(this->JXi_ + 0.5 * this->tempSum_, -xi, this->DeltaQDTH_);
            // Eq. 40
            this->delta_q_ += this->DeltaQDTH_;
            Te = getTeFK(this->Slist_, this->M_, this->delta_q_);
            xi = Log6(TransInv(Te) * this->Td_);
            // Eq. 41
            this->tempSum_.setZero();
            for (int i = 0; i < this->n_; ++i) {
                this->tempSum_ += this->DeltaQDTH_(i) * this->Hessian_[i];
            }
            getLeftInverse(this->JXi_ + this->tempSum_, -xi, this->DeltaQDMH_);
            // Eq. 42
            this->delta_q_ += this->DeltaQDMH_;
            break;

        case ALG_DTH:
            // Eq. 30
            getLeftInverse(this->JXi_, -xi, this->DeltaQDNR_);
            // Eq. 92
            getHessian(this->JXi_, this->Hessian_);
            // // Eq. 39
            this->tempSum_.setZero();
            for (int i = 0; i < this->n_; ++i) {
                this->tempSum_ += this->DeltaQDNR_(i) * this->Hessian_[i];
            }
            getLeftInverse(this->JXi_ + 0.5 * this->tempSum_, -xi, this->DeltaQDTH_);
            // Eq. 40
            this->delta_q_ += this->DeltaQDTH_;
            break;

        case ALG_DNR:
            // Eq. 30
            getLeftInverse(this->JXi_, -xi, this->DeltaQDNR_);
            this->delta_q_ += this->DeltaQDNR_;
            break;

        default:
            cout << "Invalid method specified!" << endl;
        }
    }
}




bool InverseKinematics::nearZero(double v)
{
    return std::abs(v) < EPS_NEAR_ZERO;
}

void InverseKinematics::TransToRp(const Matrix4d &T, Matrix3d &R, Vector3d &p) {
    R = T.block<3,3>(0,0);
    p = T.block<3,1>(0,3);
}

Matrix4d InverseKinematics::TransInv(const Matrix4d &T)
{
    Matrix3d R; Vector3d p;
    TransToRp(T, R, p);
    Matrix4d invT = Matrix4d::Identity();
    invT.block<3,3>(0,0) = R.transpose();
    invT.block<3,1>(0,3) = -R.transpose() * p;
    return invT;
}

// Eq. 2
Matrix3d InverseKinematics::vecToSo3(const Vector3d &v)
{
    Matrix3d m = Matrix3d::Zero();
    m(0,1) = -v(2);     m(0,2) =  v(1);
    m(1,0) =  v(2);     m(1,2) = -v(0);
    m(2,0) = -v(1);     m(2,1) =  v(0);
    return m;
}
Vector3d InverseKinematics::so3ToVec(const Matrix3d &m)
{
    Vector3d v;
    v << m(2,1), m(0,2), m(1,0);
    return v;
}

// Eq. 3
Matrix4d InverseKinematics::Exp6(const Matrix<double,6,1> &lambda)
{
    Vector3d lambda_w = lambda.segment<3>(0);
    Vector3d lambda_v = lambda.segment<3>(3);
    double n = lambda_w.norm();
    Matrix4d T = Matrix4d::Identity();
    if (nearZero(n)) {
        T.block<3,1>(0,3) = lambda_v;
        return T;
    }
    double s = sin(n/2.0) / (n/2.0);
    double c = cos(n/2.0);
    double alpha = s * c;
    double beta = s * s;
    Matrix3d lambda_w_hat = vecToSo3(lambda_w);
    Matrix3d R = Matrix3d::Identity() + alpha * lambda_w_hat + (beta/2.0) * lambda_w_hat * lambda_w_hat;
    Vector3d p = getDexpW(lambda_w) * lambda_v;
    T.block<3,3>(0,0) = R;
    T.block<3,1>(0,3) = p;
    return T;
}

//Eq. 4
Matrix3d InverseKinematics::getDexpW(const Vector3d &lambda_w)
{
    double n = lambda_w.norm();
    if (nearZero(n)) return Matrix3d::Identity();
    double s = sin(n/2.0) / (n/2.0);
    double c = cos(n/2.0);
    double alpha = s * c;
    double beta = s * s;
    Matrix3d lambda_w_hat = vecToSo3(lambda_w);
    return Matrix3d::Identity() + (beta/2.0) * lambda_w_hat + ((1.0 - alpha)/(n*n)) * (lambda_w_hat * lambda_w_hat);
}
Matrix3d InverseKinematics::getDexpWInv(const Vector3d &lambda_w) {
    double n = lambda_w.norm();
    if (nearZero(n)) return Matrix3d::Identity();
    double s = sin(n/2.0) / (n/2.0);
    double c = cos(n/2.0);
    double gamma = c / s;
    Matrix3d lambda_w_hat = vecToSo3(lambda_w);
    return Matrix3d::Identity() - 0.5 * lambda_w_hat + ((1.0 - gamma)/(n*n)) * (lambda_w_hat * lambda_w_hat);
}

// Eq. 7
Matrix<double,6,6> InverseKinematics::getDexpInv(const Matrix<double,6,1> &lambda)
{
    Vector3d lambda_w = lambda.segment<3>(0);
    Vector3d lambda_v = lambda.segment<3>(3);
    Matrix3d dexp_w_inv = getDexpWInv(lambda_w);
    double n = lambda_w.norm();
    Matrix3d lambda_w_hat = vecToSo3(lambda_w);
    Matrix3d lambda_v_hat = vecToSo3(lambda_v);
    Matrix3d D_w = Matrix3d::Zero();
    if (nearZero(n)) {
        D_w = -0.5 * lambda_v_hat;
    } else {
        double s = sin(n/2.0) / (n/2.0);
        double c = cos(n/2.0);
        double beta = s*s;
        double gamma = c / s;
        double  n2 = n * n;
        double Lambda_4 = (1.0 - gamma)/(n2);
        double Lambda_5 = (1.0/(n2)) * ((1.0/beta + gamma - 2.0)/(n2));
        D_w = -0.5 * lambda_v_hat
              + Lambda_4 * (lambda_w_hat * lambda_v_hat + lambda_v_hat * lambda_w_hat)
              + Lambda_5 * (lambda_w.dot(lambda_v)) * (lambda_w_hat * lambda_w_hat);
    }
    Matrix<double,6,6> dexp = Matrix<double,6,6>::Zero();
    dexp.block<3,3>(0,0) = dexp_w_inv;
    dexp.block<3,3>(3,0) = D_w;
    dexp.block<3,3>(3,3) = dexp_w_inv;
    return dexp;
}

// Eq. 12
Matrix<double,6,6>  InverseKinematics::adjointInv(const Matrix4d &T)
{
    Matrix3d R; Vector3d p; TransToRp(T, R, p);
    Matrix<double,6,6> AdInv = Matrix<double,6,6>::Zero();
    AdInv.block<3,3>(0,0) = R.transpose();
    AdInv.block<3,3>(0,3).setZero();
    AdInv.block<3,3>(3,0) = -R.transpose() * vecToSo3(p);
    AdInv.block<3,3>(3,3) = R.transpose();
    return AdInv;
}
Matrix<double,6,6> InverseKinematics::ad(const Matrix<double,6,1> &V)
{
    Matrix<double,6,6> A = Matrix<double,6,6>::Zero();
    Vector3d w = V.segment<3>(0);
    Vector3d v = V.segment<3>(3);
    A.block<3,3>(0,0) = vecToSo3(w);
    A.block<3,3>(3,0) = vecToSo3(v);
    A.block<3,3>(3,3) = vecToSo3(w);
    return A;
}

// Eq. 16
Matrix4d InverseKinematics::getTeFK(const MatrixXd &Slist, const MatrixXd& M, const VectorXd &delta_q)
{
    Matrix4d Te = M;
    for (int i = delta_q.size()-1; i >= 0 ; i--)
        Te = Exp6(Slist.col(i) * delta_q(i)) * Te;
    return Te;
}

// Eq. 30
void InverseKinematics::getLeftInverse(const MatrixXd &A,    const VectorXd &b,    VectorXd& x)
{
    MatrixXd ATA= A.transpose() * A;
    if (this->sigma_ > 0)
        ATA.diagonal().array() += this->sigma_;
    LLT<MatrixXd>ATA_llt(ATA);
    x =  ATA_llt.solve(A.transpose() * b);
}


// Eq. 61
Matrix<double,6,1> InverseKinematics::Log6(const Matrix4d &T)
{
    Matrix3d R; Vector3d p; TransToRp(T, R, p);
    Matrix<double,6,1> lambda = Matrix<double,6,1>::Zero();
    if ((R - Matrix3d::Identity()).norm() < EPS_NEAR_ZERO) {
        lambda.segment<3>(0).setZero();
        lambda.segment<3>(3) = p;
        return lambda;
    }
    double tr = R.trace();
    if (nearZero(tr + 1.0)) {
        // theta = pi case
        Vector3d lambda_w;
        if (!nearZero(1.0 + R(2,2))) {
            lambda_w = (PI / sqrt(2.0 * (1.0 + R(2,2)))) * Vector3d(R(0,2), R(1,2), 1.0 + R(2,2));
        } else if (!nearZero(1.0 + R(1,1))) {
            lambda_w = (PI / sqrt(2.0 * (1.0 + R(1,1)))) * Vector3d(R(0,1), 1.0 + R(1,1), R(2,1));
        } else {
            lambda_w = (PI / sqrt(2.0 * (1.0 + R(0,0)))) * Vector3d(1.0 + R(0,0), R(1,0), R(2,0));
        }
        lambda.segment<3>(0) = lambda_w;
        Matrix3d dexp_w_inv = getDexpWInv(lambda_w);
        lambda.segment<3>(3) = dexp_w_inv * p;
        return lambda;
    } else {
        double acosinput = (tr - 1.0) / 2.0;
        if (acosinput > 1.0) acosinput = 1.0;
        if (acosinput < -1.0) acosinput = -1.0;
        double theta = acos(acosinput);
        Vector3d lambda_w;
        if (nearZero(theta)) {
            lambda_w.setZero();
        } else {
            lambda_w = so3ToVec((theta/(2.0*sin(theta))) * (R - R.transpose()));
        }
        lambda.segment<3>(0) = lambda_w;
        Matrix3d dexp_w_inv = getDexpWInv(lambda_w);
        lambda.segment<3>(3) = dexp_w_inv * p;
        return lambda;
    }
}

void InverseKinematics::getTlist(const MatrixXd &Slist, const MatrixXd& M, const VectorXd &delta_q, std::vector<Matrix4d> &Tlist)
{
    Tlist[Slist.cols()] = M;
    for (int i = 0; i < Slist.cols(); ++i) {
        Tlist[i] = Exp6(Slist.col(i) * delta_q(i));
    }
}

void InverseKinematics::getTinlist(const std::vector<Matrix4d> &Tlist,  int n, std::vector<Matrix4d> &Tinlist)
{
    Tinlist[n] = Tlist[n];
    for (int i = n-1; i >= 0; --i) {
        Tinlist[i] = Tlist[i] * Tinlist[i+1];
    }
}

// Eq. 87
void InverseKinematics::getJacobian(const MatrixXd &Slist, const std::vector<Matrix4d> &Tin_list, MatrixXd &J)
{
    for (int i = 0; i < Slist.cols(); ++i)
        J.col(i) = -adjointInv(Tin_list[i+1]) * Slist.col(i);
}

// Eq. 92
void InverseKinematics::getHessian(const MatrixXd &J, vector<MatrixXd>& Hessian){
    for (int i = 1; i < J.cols(); ++i)
        Hessian[i].leftCols(i) = ad(J.col(i)) * J.leftCols(i);
}


void InverseKinematics::setAlgorithm(AlgorithmType algo){
    this->algo_ = algo;
}

void InverseKinematics::setDeltaQInit(const VectorXd &delta_q_init){
    this->delta_q_init_ = delta_q_init;
}

void InverseKinematics::setTd(Matrix4d &Td){
    this->Td_ = Td;
}

void InverseKinematics::setEpsilon(double epsilon){
    this->epsilon_ = epsilon;
}

void InverseKinematics::setTauMax(int tau_max){
    this->tau_max_ = tau_max;
}

void InverseKinematics::setSigma(double sigma){
    this->sigma_ = sigma;
}

string InverseKinematics::getAlgorithm()
{
    switch (this->algo_)
    {
        case ALG_DMH: return "ALGORITHM_DMH";
        case ALG_DTH: return "ALGORITHM_DTH";
        case ALG_DNR: return "ALGORITHM_DNR";
        default: return "UNKNOWN_ALGORITHM";
    }
}

VectorXd InverseKinematics::getDeltaQInit(){
    return this->delta_q_init_;
}

Matrix4d InverseKinematics::getTd(){
    return this->Td_;
}

double InverseKinematics::getEpsilon(){
    return this->epsilon_;
}

int InverseKinematics::getTauMax(){
    return this->tau_max_;
}

double InverseKinematics::getSigma(){
    return this->sigma_;
}

VectorXd InverseKinematics::getResult(){
    return this->delta_q_;
}

VectorXd InverseKinematics::getError(){
    return  VectorXd::Map(this->error_.data(), this->error_.size());
}


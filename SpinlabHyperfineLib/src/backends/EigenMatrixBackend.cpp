#include "pch.h"
#include "aef/backends/EigenMatrixBackend.h"

using aef::matrix::ResultCode;

// this is a really simple backend.  If only the rest were this simple 

ResultCode aef::matrix::EigenMatrixBackend::init(int argc, char** argv) {
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::shutdown() {
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::set_max_size(int nMaxDim) {
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::multiply(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
    out = A * B;
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::commutator(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
    out = A * B - B * A;
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::group_action(Eigen::MatrixXcd& out, Eigen::MatrixXcd& U, Eigen::MatrixXcd& A) {
    out = U * A * U.adjoint();
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A) {
    out = v1.adjoint() * A * v1;
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::matrix_element(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::VectorXcd& v2) {
    out = v1.adjoint() * A * v2;
    return ResultCode::Success;
}

ResultCode aef::matrix::EigenMatrixBackend::diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
    solver.compute(mat);
    evals = solver.eigenvalues();
    evecs = solver.eigenvectors();
    return ResultCode::Success;
}

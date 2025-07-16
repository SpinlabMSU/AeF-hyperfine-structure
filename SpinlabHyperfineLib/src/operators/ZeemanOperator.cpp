#include "pch.h"
#include <aef/operators/ZeemanOperator.h>

aef::operators::ZeemanOperator::ZeemanOperator(aef::MolecularSystem& sys_, double B_) : ZeemanOperator(sys_, { 0,0,B }) {
    B = B_;
}

aef::operators::ZeemanOperator::ZeemanOperator(aef::MolecularSystem& sys_, std::array<double, 3> arr) :
    sys(sys_), B_vec(arr) {
    B = sqrt(B_x * B_x + B_y * B_y + B_z * B_z);
    info.name = fmt::format("Zeeman_({},{},{})", B_x, B_y, B_z);
    info.description = fmt::format("Zeeman Shift, electric field is ({}, {}, {}) V/cm", B_x, B_y, B_z);
    info.is_hermitian = 1;
    info.is_potential = 1;
}



aef::operators::ZeemanOperator::~ZeemanOperator() {
}

dcomplex aef::operators::ZeemanOperator::matrixElement(size_t kdx1, size_t kdx2) {
    std::array<dcomplex, 3> mdm = sys.get_calc()->molec_mdm(kdx1, kdx2);
    dcomplex melt = 0;
    for (int idx = 0; idx < 3; idx++) {
        melt += B_vec[idx] * mdm[idx];
    }
    return melt;
}

void aef::operators::ZeemanOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {
    matrix.setZero();
    constexpr double eps = 1e-8;
    constexpr double eps_sq = eps * eps;
    for (int jdx = 0; jdx < matrix.cols(); jdx++) {
        for (int idx = 0; idx < matrix.rows(); idx++) {
            dcomplex melt = matrixElement(idx, jdx);
            // only insert significantly nonzero elements to aid in preserving sparsity
            if (std::norm(melt) > eps_sq) {
                matrix.coeffRef(idx, jdx) = matrixElement(idx, jdx);
            }
        }
    }
}

void aef::operators::ZeemanOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    matrix.setZero();
    for (int jdx = 0; jdx < matrix.cols(); jdx++) {
        for (int idx = 0; idx < matrix.rows(); idx++) {
            matrix(idx, jdx) = matrixElement(idx, jdx);
        }
    }
}

aef::operators::OperatorInfo* aef::operators::ZeemanOperator::getInfo() {
    return &info;
}

// ket version

aef::operators::ket::ZeemanOperator::ZeemanOperator(double B_) : ZeemanOperator({0,0,B}) {
    B = B_;
}

aef::operators::ket::ZeemanOperator::ZeemanOperator(std::array<double, 3> arr) :
    B_vec(arr)
{
    B = sqrt(B_x * B_x + B_y * B_y + B_z * B_z);
    info.name = fmt::format("Zeeman_({},{},{})", B_x, B_y, B_z);
    info.description = fmt::format("Zeeman Shift, electric field is ({}, {}, {}) V/cm", B_x, B_y, B_z);
    info.is_hermitian = 1;
    info.is_potential = 1;
}



aef::operators::ket::ZeemanOperator::~ZeemanOperator() {}

dcomplex aef::operators::ket::ZeemanOperator::matrixElement(basis_ket k1, basis_ket k2) {
    return dcomplex();
}

void aef::operators::ket::ZeemanOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {

    constexpr double eps = 1e-8;
    constexpr double eps_sq = eps * eps;

    Eigen::Index num_significant = 0;

    // j first because column-major 
    for (Eigen::Index j = 0; j < matrix.cols(); j++) {
        basis_ket kj = basis_ket::from_index(j);
        for (Eigen::Index i = 0; i < matrix.rows(); i++) {
            basis_ket ki = basis_ket::from_index(i);
            std::array<dcomplex, 3> mdm = ki.molec_mdm(kj);
            // V_stk = - \vec{\mu} \cdot \vec{E}
            dcomplex mat_elt = -(B_x * mdm[0] + B_y * mdm[1] + B_z * mdm[2]);

            if (std::norm(mat_elt) > eps_sq) {
                // fill matrix
                matrix.coeffRef(i, j) = mat_elt;
                num_significant++;
            }

        }
    }
}

void aef::operators::ket::ZeemanOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    // j first because column-major 
    for (Eigen::Index j = 0; j < matrix.cols(); j++) {
        basis_ket kj = basis_ket::from_index(j);
        for (Eigen::Index i = 0; i < matrix.rows(); i++) {
            basis_ket ki = basis_ket::from_index(i);
            std::array<dcomplex, 3> mdm = ki.molec_mdm(kj);
            // V_stk = - \vec{\mu} \cdot \vec{E}
            matrix(i, j) = -(B_x * mdm[0] + B_y * mdm[1] + B_z * mdm[2]);
        }
    }
}

aef::operators::OperatorInfo* aef::operators::ket::ZeemanOperator::getInfo() {
    return &info;
}

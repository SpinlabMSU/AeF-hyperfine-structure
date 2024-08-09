#include "pch.h"
#include <aef/operators/StarkOperator.h>

aef::operators::StarkOperator::StarkOperator(double E_, Direction dir_)
    : E(E_), dir(dir_), E_vec({ 0,0,0 }) {
    if (dir < Direction::X || dir > Direction::Z) {
        dir = Direction::X;
    }
    int idir = static_cast<int>(dir);
    info.name = fmt::format("Stark_{}", dirString(dir));
    info.description = fmt::format("Stark Shift field along {} axis", dirString(dir));
    info.is_hermitian = 1;
    info.is_potential = 1;
    E_vec[idir] = E_;
}

aef::operators::StarkOperator::StarkOperator(std::array<double, 3> E_vec_): 
    E_vec(E_vec_), dir(Direction::INVALID) {
    E = std::sqrt(E_x * E_x + E_y * E_y + E_z * E_z);
    info.name = fmt::format("Stark_({},{},{})", E_x, E_y, E_z);
    info.description = fmt::format("Stark Shift, electric field is ({}, {}, {}) V/cm", E_x, E_y, E_z);
    info.is_hermitian = 1;
    info.is_potential = 1;
}

aef::operators::StarkOperator::~StarkOperator() {
}

dcomplex aef::operators::StarkOperator::matrixElement(basis_ket ki, basis_ket kj) {
    std::array<dcomplex, 3> stk = ki.molec_edm(kj);
    return E_x * stk[0] + E_y * stk[1] + E_z * stk[2];
}

void aef::operators::StarkOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {
}

void aef::operators::StarkOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    // j first because 
    for (int j = 0; j < matrix.cols(); j++) {
        basis_ket kj = basis_ket::from_index(j);
        for (int i = 0; i < matrix.rows(); i++) {
            basis_ket ki = basis_ket::from_index(i);
            std::array<dcomplex, 3> stk = ki.molec_edm(kj);
            matrix(i, j) = E_x * stk[0] + E_y * stk[1] + E_z * stk[2];
        }
    }
}

aef::operators::OperatorInfo* aef::operators::StarkOperator::getInfo() {
    return &info;
}

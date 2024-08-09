#include "pch.h"
#include <aef/operators/StarkOperator.h>

aef::operators::StarkOperator::StarkOperator(double E_, Direction dir_)
: E(E_), dir(dir_) {
    info.name = fmt::format("Stark_{}", dirString(dir_));
    info.description = fmt::format("Stark Shift field along {} axis", dirString(dir));
    info.is_hermitian = 1;
    info.is_potential = 1;
}

aef::operators::StarkOperator::~StarkOperator() {
}

dcomplex aef::operators::StarkOperator::matrixElement(basis_ket k1, basis_ket k2) {
    return dcomplex();
}

void aef::operators::StarkOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {
}

void aef::operators::StarkOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    // j first because 
    for (int j = 0; j < matrix.cols(); j++) {
        basis_ket kj = basis_ket::from_index(j);
        for (int i = 0; i < matrix.rows(); i++) {
            basis_ket ki = basis_ket::from_index(i);
            dcomplex elt = 0;
            using namespace std::complex_literals;
            constexpr double inv_sqrt2 = std::numbers::sqrt2 / 2.0;
            if (dir == Direction::X) {
                elt = (ki.d1t(kj) - ki.d11(kj)) * inv_sqrt2;
            } else if (dir == Direction::Y) {
                elt = (ki.d1t(kj) + ki.d11(kj)) * inv_sqrt2 * 1i;
            } else {
                elt = ki.d10(kj);
            }
            constexpr double mu = 1;// TODO fetch this from ki/kj
            matrix(i, j) = elt * E * mu;
        }
    }
}

aef::operators::OperatorInfo* aef::operators::StarkOperator::getInfo() {
    return &info;
}

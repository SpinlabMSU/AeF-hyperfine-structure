#include "pch.h"
#include <aef/operators/eEDMOperator.h>

aef::operators::eEDMOperator::eEDMOperator() {
    // TODO
}

aef::operators::eEDMOperator::~eEDMOperator() {}

dcomplex aef::operators::eEDMOperator::matrixElement(basis_ket k1, basis_ket k2) {
    return dcomplex();
}

void aef::operators::eEDMOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {}

void aef::operators::eEDMOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    // j first because 
    for (int j = 0; j < matrix.cols(); j++) {
        basis_ket kj = basis_ket::from_index(j);
        for (int i = 0; i < matrix.rows(); i++) {
            basis_ket ki = basis_ket::from_index(i);
            dcomplex elt = 0;
            using namespace std::complex_literals;
            constexpr double inv_sqrt2 = std::numbers::sqrt2 / 2.0;
            matrix(i, j) = 0;//TODO implement
        }
    }
}

aef::operators::OperatorInfo *aef::operators::eEDMOperator::getInfo() {
    return &info;
}

#include "pch.h"
#include <aef/operators/NSMOperator.h>

aef::operators::NSMOperator::NSMOperator() {
    info.name = std::string("NSM");
    info.description = std::string("Nuclear Schiff Moment Operator");
    info.is_hermitian = 1;
    info.is_potential = 1;

    //
}

aef::operators::NSMOperator::~NSMOperator() {
    memset(((void*)this), 0, sizeof(*this));
}

dcomplex aef::operators::NSMOperator::matrixElement(basis_ket k1, basis_ket k2) {
    return 0;
}

void aef::operators::NSMOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {
}

void aef::operators::NSMOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    // j first because 
    for (int j = 0; j < matrix.cols(); j++) {
        basis_ket kj = basis_ket::from_index(j);
        for (int i = 0; i < matrix.rows(); i++) {
            basis_ket ki = basis_ket::from_index(i);
            // TODO implement NSM on light nucleus
            matrix(i, j) = 0;
        }
    }
}

aef::operators::OperatorInfo *aef::operators::NSMOperator::getInfo() {
    return &info;
}

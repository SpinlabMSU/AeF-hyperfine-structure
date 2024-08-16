#include "pch.h"
#include <aef/operators/ZeemanOperator.h>

aef::operators::ZeemanOperator::ZeemanOperator() {}

aef::operators::ZeemanOperator::~ZeemanOperator() {}

dcomplex aef::operators::ZeemanOperator::matrixElement(basis_ket k1, basis_ket k2) {
    return dcomplex();
}

void aef::operators::ZeemanOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {}

void aef::operators::ZeemanOperator::fillMatrix(Eigen::MatrixXcd& matrix) {}

aef::operators::OperatorInfo* aef::operators::ZeemanOperator::getInfo() {
    return nullptr;
}

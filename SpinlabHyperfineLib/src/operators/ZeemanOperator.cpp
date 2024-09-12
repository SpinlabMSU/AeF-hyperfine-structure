#include "pch.h"
#include <aef/operators/ZeemanOperator.h>

aef::operators::ZeemanOperator::ZeemanOperator(double B_) : ZeemanOperator({0,0,B})
{
    B = B_;
}

aef::operators::ZeemanOperator::ZeemanOperator(std::array<double, 3> arr) :
    B_vec(arr)
{
    B = sqrt(B_x * B_x + B_y * B_y + B_z * B_z);
    info.name = fmt::format("Zeeman_({},{},{})", B_x, B_y, B_z);
    info.description = fmt::format("Zeeman Shift, electric field is ({}, {}, {}) V/cm", B_x, B_y, B_z);
    info.is_hermitian = 1;
    info.is_potential = 1;
}



aef::operators::ZeemanOperator::~ZeemanOperator() {}

dcomplex aef::operators::ZeemanOperator::matrixElement(basis_ket k1, basis_ket k2) {
    return dcomplex();
}

void aef::operators::ZeemanOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {}

void aef::operators::ZeemanOperator::fillMatrix(Eigen::MatrixXcd& matrix) {}

aef::operators::OperatorInfo* aef::operators::ZeemanOperator::getInfo() {
    return &info;
}

#include "pch.h"
#include <aef/operators/NSMOperator.h>

aef::operators::NSMOperator::NSMOperator(aef::MolecularSystem& sys_, bool on_heavy_nucleus_): sys(sys_), on_heavy_nucleus(on_heavy_nucleus_) {
    if (on_heavy_nucleus) {
        info.name = std::string("NSM_heavy_angdep");
        info.description = std::string("Heavy-nucleus Nuclear Schiff Moment Angular-Dependence Operator");
    } else {
        info.name = std::string("NSM_light_angdep");
        info.description = std::string("Light-nucleus Nuclear Schiff Moment Angular-Dependence Operator");
    }
    info.is_hermitian = 1;
    info.is_potential = 1;
}

aef::operators::NSMOperator::~NSMOperator() {
}

dcomplex aef::operators::NSMOperator::matrixElement(size_t kdx1, size_t kdx2) {
    return nan;
}

void aef::operators::NSMOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {
}

void aef::operators::NSMOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    if (on_heavy_nucleus) {
        sys.get_calc()->calculate_I1_dot_ina(matrix);
    } else {
        sys.get_calc()->calculate_I2_dot_ina(matrix);
    }
}

aef::operators::OperatorInfo* aef::operators::NSMOperator::getInfo() {
    return &info;
}




aef::operators::ket::NSMOperator::NSMOperator() {
    info.name = std::string("NSM_angdep");
    info.description = std::string("Light-nucleus Nuclear Schiff Moment Angular-Dependence Operator");
    info.is_hermitian = 1;
    info.is_potential = 1;

    //
}

aef::operators::ket::NSMOperator::~NSMOperator() {
}

dcomplex aef::operators::ket::NSMOperator::matrixElement(basis_ket k1, basis_ket k2) {
    return k1.I_dot_ina(k2);
}

void aef::operators::ket::NSMOperator::fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix) {
}

void aef::operators::ket::NSMOperator::fillMatrix(Eigen::MatrixXcd& matrix) {
    // j first because 
    for (int j = 0; j < matrix.cols(); j++) {
        basis_ket kj = basis_ket::from_index(j);
        for (int i = 0; i < matrix.rows(); i++) {
            basis_ket ki = basis_ket::from_index(i);
            matrix(i, j) = ki.I_dot_ina(kj);
        }
    }
}

aef::operators::OperatorInfo *aef::operators::ket::NSMOperator::getInfo() {
    return &info;
}

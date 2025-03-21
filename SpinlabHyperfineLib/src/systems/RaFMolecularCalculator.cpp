#include "pch.h"
#include "aef/systems/RaFMolecularCalculator.h"
using aef::ResultCode;


aef::RaFMolecularCalculator::RaFMolecularCalculator(): nmax(-1), nBasisElts(0){
}

aef::RaFMolecularCalculator::RaFMolecularCalculator(spin nmax_) : nmax(nmax_) {
    set_nmax(nmax_);
}

aef::RaFMolecularCalculator::~RaFMolecularCalculator() {}

ResultCode aef::RaFMolecularCalculator::get_parameter(std::string id, double& out) {
    return ResultCode();
}

ResultCode aef::RaFMolecularCalculator::set_parameter(std::string id, double value) {
    return ResultCode();
}

spin aef::RaFMolecularCalculator::get_nmax() {
    return nmax;
}

void aef::RaFMolecularCalculator::set_nmax(spin nmax_) {
    nmax = nmax_;
    nBasisElts = jf_basis_vec::index_of_n(nmax_ + 1);
    // construct basis
    basis.clear();
    basis.reserve(nBasisElts);
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        basis.emplace_back(jf_basis_vec::from_index(idx));
    }
}

aef::universal_diatomic_basis_vec aef::RaFMolecularCalculator::get_basis_ket(int idx) {
    if (idx >= nBasisElts) {
        double nan = std::nan("5555");
        return universal_diatomic_basis_vec(nan, nan, nan, nan, nan);
    }
    jf_basis_vec v = basis[idx];
    return universal_diatomic_basis_vec(v.n, v.j, v.f1, v.f, v.m_f);
}

int aef::RaFMolecularCalculator::get_index(universal_diatomic_basis_vec v) {
    jf_basis_vec ket(v.n, v.j, v.f_1, v.f, v.m_f);
    return ket.index();
}

ResultCode aef::RaFMolecularCalculator::calculate_H_rot(Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic>& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        H.diagonal()(idx, idx) = basis[idx].H_rot();
    }
    return ResultCode::Success;
}

ResultCode aef::RaFMolecularCalculator::calculate_H_hfs(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].H_hfs(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
    return ResultCode::Success;
}

ResultCode aef::RaFMolecularCalculator::calculate_H_dev(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].H_hfs(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
    return ResultCode::Success;;
}

ResultCode aef::RaFMolecularCalculator::calculate_H_stk(Eigen::MatrixXcd& H, double E_z) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = E_z * basis[idx].H_hfs(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
    return ResultCode::Success;;
}

void aef::RaFMolecularCalculator::calculate_F_z(Eigen::MatrixXcd& F_z) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        F_z(idx, idx) = basis[idx].m_f;
    }
}

ResultCode aef::RaFMolecularCalculator::calculate_dkq(Eigen::MatrixXcd& d, int q) {
    if (q == -1) {
        calculate_d1t(d);
        return ResultCode::Success;
    } else  if (q == 0) {
        calculate_d10(d);
        return ResultCode::Success;
    } else if (q == 1) {
        calculate_d11(d);
        return ResultCode::Success;
    }
    return ResultCode::InvalidArgument;
}

void aef::RaFMolecularCalculator::calculate_d1t(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].d1t(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::RaFMolecularCalculator::calculate_d10(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].d10(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::RaFMolecularCalculator::calculate_d11(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].d11(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::RaFMolecularCalculator::load(std::istream& in) {
    return;
}

void aef::RaFMolecularCalculator::save(std::ostream& out) {
    return;
}

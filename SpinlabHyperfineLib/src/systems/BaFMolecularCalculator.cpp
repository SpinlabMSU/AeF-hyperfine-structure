#include "pch.h"
#include "aef/systems/BaFMolecularCalculator.h"
using aef::ResultCode;


aef::BaFMolecularCalculator::BaFMolecularCalculator(): nmax(-1), nBasisElts(0){
}

aef::BaFMolecularCalculator::BaFMolecularCalculator(spin nmax_) : nmax(nmax_) {
    set_nmax(nmax_);
}

aef::BaFMolecularCalculator::~BaFMolecularCalculator() {}

ResultCode aef::BaFMolecularCalculator::get_parameter(std::string id, double& out) {
    return ResultCode::Unimplemented;
}

ResultCode aef::BaFMolecularCalculator::set_parameter(std::string id, double value) {
    return ResultCode::Unimplemented;
}

spin aef::BaFMolecularCalculator::get_nmax() {
    return nmax;
}

void aef::BaFMolecularCalculator::set_nmax(spin nmax_) {
    nmax = nmax_;
    nBasisElts = j_basis_vec::index_of_n(nmax_ + 1);
    // construct basis
    basis.clear();
    basis.reserve(nBasisElts);
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        basis.emplace_back(j_basis_vec::from_index(idx));
    }
    using aef::half;
    lowest_states = {
        aef::universal_diatomic_basis_vec(0, half, 0,  0),
        aef::universal_diatomic_basis_vec(0, half, 1, -1),
        aef::universal_diatomic_basis_vec(0, half, 1,  0),
        aef::universal_diatomic_basis_vec(0, half, 1, +1),
    };
}

#include <numeric>

aef::universal_diatomic_basis_vec aef::BaFMolecularCalculator::get_basis_ket(int idx) {
    if (idx >= nBasisElts) {
        double nan = std::nan("5555");
        return universal_diatomic_basis_vec(nan, nan, nan, nan, nan);
    }
    j_basis_vec v = basis[idx];
    return universal_diatomic_basis_vec(v.n, v.j, v.f, v.m_f);
}

int aef::BaFMolecularCalculator::get_index(universal_diatomic_basis_vec v) {
    j_basis_vec ket(v.n, v.j, v.f, v.m_f);
    return ket.index();
}

ResultCode aef::BaFMolecularCalculator::calculate_H_rot(Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic>& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        H.diagonal()[idx] = basis[idx].H_rot();
    }
    return ResultCode::Success;
}

ResultCode aef::BaFMolecularCalculator::calculate_H_hfs(Eigen::MatrixXcd& H) {
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

ResultCode aef::BaFMolecularCalculator::calculate_H_dev(Eigen::MatrixXcd& H, double K) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].H_dev(basis[jdx], K);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
    return ResultCode::Success;;
}

ResultCode aef::BaFMolecularCalculator::calculate_H_stk(Eigen::MatrixXcd& H, double E_z) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].H_st(basis[jdx], E_z);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
    return ResultCode::Success;;
}

void aef::BaFMolecularCalculator::calculate_F_z(Eigen::MatrixXcd& F_z) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        F_z(idx, idx) = basis[idx].m_f;
    }
}

ResultCode aef::BaFMolecularCalculator::calculate_dkq(Eigen::MatrixXcd& d, int q) {
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

void aef::BaFMolecularCalculator::calculate_d1t(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].d1t(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::BaFMolecularCalculator::calculate_d10(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].d10(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::BaFMolecularCalculator::calculate_d11(Eigen::MatrixXcd& H) {
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].d11(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::BaFMolecularCalculator::load(std::istream& in) {
    return;
}

void aef::BaFMolecularCalculator::save(std::ostream& out) {
    return;
}

const char* aef::BaFMolecularCalculator::get_calc_type() {
    return calc_type_str;
}

int aef::BaFMolecularCalculator::get_lowest_rotational_state_size() {
    return aef::j_basis_vec::index_of_n(1);
}

std::vector<universal_diatomic_basis_vec> aef::BaFMolecularCalculator::get_lowest_states() {
    return lowest_states;
}

std::array<dcomplex, 3> aef::BaFMolecularCalculator::molec_edm(int kdx1, int kdx2) {
    return basis[kdx1].molec_edm(basis[kdx2]);
}

std::array<dcomplex, 3> aef::BaFMolecularCalculator::molec_mdm(int kdx1, int kdx2) {
    return basis[kdx1].molec_mdm(basis[kdx2]);
}

void aef::BaFMolecularCalculator::calculate_S_dot_ina(Eigen::MatrixXcd& A) {
    A.setZero();
    for (int jdx = 0; jdx < nBasisElts; jdx++) {
        for (int idx = 0; idx <= jdx; idx++) {
            dcomplex melt = basis[idx].S_dot_ina(basis[jdx]);
            A(idx, jdx) = melt;
            A(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::BaFMolecularCalculator::calculate_I1_dot_ina(Eigen::MatrixXcd& A) {
    // I1 corresponds to the heavy nucleus, which has spin zero
    A.setZero();
}

void aef::BaFMolecularCalculator::calculate_I2_dot_ina(Eigen::MatrixXcd& A) {
    // the light nucleus is associated with spin I_2, technically.
    A.setZero();
    for (int jdx = 0; jdx < nBasisElts; jdx++) {
        for (int idx = 0; idx <= jdx; idx++) {
            dcomplex melt = basis[idx].I_dot_ina(basis[jdx]);
            A(idx, jdx) = melt;
            A(jdx, idx) = std::conj(melt);
        }
    }
}
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
    return ResultCode::Unimplemented;
}

ResultCode aef::RaFMolecularCalculator::set_parameter(std::string id, double value) {
    return ResultCode::Unimplemented;
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
    std::cout << "Idx, basis ket" << std::endl;
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        basis.emplace_back(jf_basis_vec::from_index(idx));
        std::cout << idx << ", " << basis[idx].ket_csv_str() << std::endl;
    }
    std::cout << std::endl << std::endl;
    using aef::half;
    lowest_states = {
        // f1 = 0, f = 1/2 doublet
        aef::universal_diatomic_basis_vec(0, half, 0, half, -half),
        aef::universal_diatomic_basis_vec(0, half, 0, half, +half),

        // f1 = 1, f = 1/2 doublet
        aef::universal_diatomic_basis_vec(0, half, 1, 1/2., -1/2.),
        aef::universal_diatomic_basis_vec(0, half, 1, 1/2., +1/2.),

        // f1 = 1, f = 3/2 quadruplet
        aef::universal_diatomic_basis_vec(0, half, 1, 3/2., -3 / 2.),
        aef::universal_diatomic_basis_vec(0, half, 1, 3/2., -1 / 2.),
        aef::universal_diatomic_basis_vec(0, half, 1, 3/2., +1 / 2.),
        aef::universal_diatomic_basis_vec(0, half, 1, 3/2., +3 / 2.),
    };
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
    H.setZero();
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        H.diagonal()[idx] = basis[idx].H_rot();
    }
    return ResultCode::Success;
}

ResultCode aef::RaFMolecularCalculator::calculate_H_hfs(Eigen::MatrixXcd& H) {
    H.setZero();
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

ResultCode aef::RaFMolecularCalculator::calculate_H_dev(Eigen::MatrixXcd& H, double K) {
    H.setZero();
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

ResultCode aef::RaFMolecularCalculator::calculate_H_stk(Eigen::MatrixXcd& H, double E_z) {
    H.setZero();
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // optimization removed temporarily -- prev: operators are hermitian matricies
        for (size_t jdx = 0; jdx < nBasisElts; jdx++) {
            dcomplex melt = basis[idx].H_st(basis[jdx], E_z);
            H(idx, jdx) = melt;
            //H(jdx, idx) = std::conj(melt);
        }
    }
    return ResultCode::Success;;
}

void aef::RaFMolecularCalculator::calculate_F_z(Eigen::MatrixXcd& F_z) {
    F_z.setZero();
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        F_z(idx, idx) = basis[idx].m_f;
    }
}

ResultCode aef::RaFMolecularCalculator::calculate_dkq(Eigen::MatrixXcd& d, int q) {
    d.setZero();
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
    H.setZero();
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
    H.setZero();
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
    H.setZero();
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        // operators are hermitian matricies
        for (size_t jdx = 0; jdx <= idx; jdx++) {
            dcomplex melt = basis[idx].d11(basis[jdx]);
            H(idx, jdx) = melt;
            H(jdx, idx) = std::conj(melt);
        }
    }
}

std::array<dcomplex, 3> aef::RaFMolecularCalculator::molec_edm(int kdx1, int kdx2) {
    return basis[kdx1].molec_edm(basis[kdx2]);
}

std::array<dcomplex, 3> aef::RaFMolecularCalculator::molec_mdm(int kdx1, int kdx2) {
    return basis[kdx1].molec_mdm(basis[kdx2]);
}

void aef::RaFMolecularCalculator::calculate_S_dot_ina(Eigen::MatrixXcd& A) {
    A.setZero();
    for (int jdx = 0; jdx < nBasisElts; jdx++) {
        for (int idx = 0; idx <= jdx; idx++) {
            dcomplex melt = basis[idx].S_dot_ina(basis[jdx]);
            A(idx, jdx) = melt;
            A(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::RaFMolecularCalculator::calculate_I1_dot_ina(Eigen::MatrixXcd& A) {
    A.setZero();
    for (int jdx = 0; jdx < nBasisElts; jdx++) {
        for (int idx = 0; idx <= jdx; idx++) {
            dcomplex melt = basis[idx].I1_dot_ina(basis[jdx]);
            A(idx, jdx) = melt;
            A(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::RaFMolecularCalculator::calculate_I2_dot_ina(Eigen::MatrixXcd& A) {
    A.setZero();
    for (int jdx = 0; jdx < nBasisElts; jdx++) {
        for (int idx = 0; idx <= jdx; idx++) {
            dcomplex melt = basis[idx].I2_dot_ina(basis[jdx]);
            A(idx, jdx) = melt;
            A(jdx, idx) = std::conj(melt);
        }
    }
}

void aef::RaFMolecularCalculator::load(std::istream& in) {
    return;
}

void aef::RaFMolecularCalculator::save(std::ostream& out) {
    return;
}

const char* aef::RaFMolecularCalculator::get_calc_type() {
    return calc_type_str;
}

int aef::RaFMolecularCalculator::get_lowest_rotational_state_size() {
    return aef::jf_basis_vec::index_of_n(1);
}

std::vector<universal_diatomic_basis_vec> aef::RaFMolecularCalculator::get_lowest_states() {
    return lowest_states;
}

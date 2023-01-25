#include "HyperfineCalculator.hpp"

HyperfineCalculator::HyperfineCalculator(spin nmax_, double E_z_, bool enableDev_)
    : nmax(nmax_), E_z(E_z_), init(false), enableDev(enableDev_), nBasisElts(j_basis_vec::index_of_n(nmax_ + 1)), diagonalized(false) {

    // construct basis
    basis.reserve(nBasisElts);
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        basis.emplace_back(j_basis_vec::from_index(idx));
    }

    //this->H_rot = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    H_rot.resize(nBasisElts);
    this->H_hfs = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    this->H_stk = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    this->H_dev = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    this->H_tot = Eigen::MatrixXcd(nBasisElts, nBasisElts);
}

HyperfineCalculator::~HyperfineCalculator() {
}

bool HyperfineCalculator::calculate_matrix_elts() {

    // rotational hamiltonian
    for (int idx = 0; idx < nBasisElts; idx++) {
        H_rot.diagonal()[idx] = basis[idx].H_rot();
    }

    // stark shift
    for (int idx = 0; idx < nBasisElts; idx++) {
        for (int jdx = 0; jdx < nBasisElts; jdx++) {
            H_stk(idx, jdx) = basis[idx].H_st(basis[jdx], E_z);
        }
    }

    // hyperfine shift:
    for (int idx = 0; idx < nBasisElts; idx++) {
        for (int jdx = 0; jdx < nBasisElts; jdx++) {
            H_hfs(idx, jdx) = basis[idx].H_hfs(basis[jdx]);
        }
    }

    // devonshire
    if (enableDev) {
        for (int idx = 0; idx < nBasisElts; idx++) {
            for (int jdx = 0; jdx < nBasisElts; jdx++) {
                H_dev(idx, jdx) = basis[idx].H_hfs(basis[jdx]);
            }
        }
    } else {
        H_dev.setZero();
    }

    H_tot.setZero();
    H_tot.diagonal() = H_rot.diagonal();
    H_tot += H_stk + H_hfs + H_dev;

    init = true;
    diagonalized = false;
    return init;
}

bool HyperfineCalculator::diagonalize_H() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
    solver.compute(H_tot);
    Es = solver.eigenvalues();
    Vs = solver.eigenvectors();
    diagonalized = true;
    return diagonalized;
}

bool HyperfineCalculator::load_matrix_elts(std::istream& in) {
    return false;
}

bool HyperfineCalculator::save_matrix_elts(std::ostream& out) {
    return false;
}

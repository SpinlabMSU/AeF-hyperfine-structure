#include "HyperfineCalculator.hpp"

HyperfineCalculator::HyperfineCalculator(spin nmax_, double E_z_, bool enableDev_)
    : nmax(nmax_), E_z(E_z_), init(false), enableDev(enableDev_), nBasisElts(j_basis_vec::index_of_n(nmax_ + 1)), diagonalized(false) {

    set_nmax(nmax_);

    //this->H_rot = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    H_rot.resize(nBasisElts);
    this->H_hfs = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    this->H_stk = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    this->H_dev = Eigen::MatrixXcd(nBasisElts, nBasisElts);
    this->H_tot = Eigen::MatrixXcd(nBasisElts, nBasisElts);
}

void HyperfineCalculator::set_nmax(spin nmax_) {
    nmax = nmax_;
    nBasisElts = j_basis_vec::index_of_n(nmax_ + 1);
    // construct basis
    basis.reserve(nBasisElts);
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        basis.emplace_back(j_basis_vec::from_index(idx));
    }
    init = false;
    diagonalized = false;

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

#include "matrix_io.hpp"
#include <zstr.hpp>
#include <algorithm>
#include <iterator>

static constexpr uint8_t MAGIC[8] = { 'A', 'e', 'F', 0, 'H', 'D', 'a', 't' };
static constexpr uint16_t VERSION = 0;

enum hyp_flags : uint16_t {
    FLAG_DIAG = 1
};


bool HyperfineCalculator::load_matrix_elts(std::istream& in) {
    zstr::istream zin(in);
    //auto& zin = in;
    uint8_t test[8] = {};

    std::cout << "READING MAGIC :";
    for (int idx = 0; idx < sizeof(test); idx++) {
        zin >> test[idx];
        std::cout << test[idx];
    }
    std::cout << std::endl;
    //zin.read((char*)test, 8);
    if (memcmp(test, MAGIC, sizeof(MAGIC)) != 0) {
        union {
            uint8_t* u8;
            uint64_t* u64;
        } data;
        data.u8 = test;
        uint64_t rmag = *data.u64;
        data.u8 = (uint8_t*)MAGIC;
        std::cout << "BAD MAGIC " << rmag << " GOOD WOULD BE " << *data.u64 << std::endl;
        return false;
    } else {
        union {
            uint8_t* u8;
            uint64_t* u64;
        } data;
        data.u8 = test;
        std::cout << "READ CORRECT MAGIC " << *data.u64 << std::endl;
    }
    uint16_t version = 0;
    zin.read((char*)&version, sizeof(version));

    if (version > VERSION) {
        std::cout << "BAD VERSION " << version << std::endl;
        uint16_t flags;
        zin.read((char*)&flags, sizeof(flags));
        std::cout << "FLAGS WOULD BE " << flags << std::endl;
        return false;
    }

    uint16_t flags;
    zin.read((char*)&flags, sizeof(flags));

    
    diagonalized = flags & FLAG_DIAG;

    uint32_t nmax_ = 0;
    zin.read((char*)&nmax_, sizeof(nmax_));
    std::cout << "NMax is " << nmax_ << std::endl;
    nmax = nmax_;
    set_nmax(nmax_);
    
    Eigen::read_binary(zin, H_rot.diagonal());
    Eigen::read_binary(zin, H_hfs);
    Eigen::read_binary(zin, H_stk);
    Eigen::read_binary(zin, H_dev);
    Eigen::read_binary(zin, H_tot);
    
    init = true;

    if (diagonalized) {
        Eigen::read_binary(zin, Es);
        Eigen::read_binary(zin, Vs);
    }

    return true;
}

bool HyperfineCalculator::save_matrix_elts(std::ostream& out) {
    zstr::ostream zout(out);
    //auto& zout = out;
    std::copy(MAGIC, MAGIC + 8, std::ostream_iterator<uint8_t>(zout));
    zout.write((char*)&VERSION, sizeof(VERSION));
    uint16_t flags = 0;

    if (diagonalized) {
        flags |= FLAG_DIAG;
    }

    zout.write((char*)&flags, sizeof(flags));
    uint32_t nmax_ = (uint32_t)nmax;
    zout.write((char*)&nmax_, sizeof(nmax_));

    Eigen::write_binary(zout, H_rot.diagonal());
    Eigen::write_binary(zout, H_hfs);
    Eigen::write_binary(zout, H_stk);
    Eigen::write_binary(zout, H_dev);
    Eigen::write_binary(zout, H_tot);

    if (diagonalized) {
        Eigen::write_binary(zout, Es);
        Eigen::write_binary(zout, Vs);
    }
    
    return true;
}
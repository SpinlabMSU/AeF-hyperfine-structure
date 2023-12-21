/*
    This file is part of the AeF-hyperfine-structure program. 
    
    AeF-hyperfine-structure is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version.

    AeF-hyperfine-structure is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    AeF-hyperfine-structure. If not, see <https://www.gnu.org/licenses/>.
*/
#include <pch.h>
#include "aef/aef.h"
#include <bitset>

namespace fs = std::filesystem;

HyperfineCalculator::HyperfineCalculator(spin nmax_, double E_z_, double K_)
    : nmax(nmax_), E_z(E_z_), init(false), enableDev(K_ != 0),
      nBasisElts(j_basis_vec::index_of_n(nmax_ + 1)), diagonalized(false),
      K(K_), dkq_init(false), load_version(HyperfineCalculator::CURRENT_VERSION){

  set_nmax(nmax_);
}

void HyperfineCalculator::set_nmax(spin nmax_) {
  nmax = nmax_;
  nBasisElts = j_basis_vec::index_of_n(nmax_ + 1);
  // construct basis
  basis.clear();
  basis.reserve(nBasisElts);
  for (size_t idx = 0; idx < nBasisElts; idx++) {
    basis.emplace_back(j_basis_vec::from_index(idx));
  }

  for (int idx = 0; idx < nBasisElts; idx++) {
    assert(basis[idx] == j_basis_vec::from_index(idx));
  }
  init = false;
  diagonalized = false;

  H_rot.resize(nBasisElts);
  this->H_hfs = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  this->H_stk = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  this->H_dev = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  this->H_tot = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  this->F_z = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  this->d10 = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  this->d11 = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  this->d1t = Eigen::MatrixXcd(nBasisElts, nBasisElts);
  Es.resize(nBasisElts);
  Vs.resize(nBasisElts, nBasisElts);
  Vst.resize(nBasisElts, nBasisElts);
}

HyperfineCalculator::~HyperfineCalculator() {}

bool HyperfineCalculator::calculate_matrix_elts() {

  // rotational hamiltonian and F_z
  for (int idx = 0; idx < nBasisElts; idx++) {
    H_rot.diagonal()[idx] = basis[idx].H_rot();
    F_z(idx, idx) = basis[idx].m_f;
  }

  // stark shift
  for (int idx = 0; idx < nBasisElts; idx++) {
    for (int jdx = 0; jdx < nBasisElts; jdx++) {
      H_stk(idx, jdx) = basis[idx].H_st(basis[jdx], E_z);
      d10(idx, jdx) = basis[idx].d10(basis[jdx]);
      d11(idx, jdx) = basis[idx].d11(basis[jdx]);
      d1t(idx, jdx) = basis[idx].d1t(basis[jdx]);
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
        H_dev(idx, jdx) = basis[idx].H_dev(basis[jdx], K);
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

bool HyperfineCalculator::diagonalize_H(bool use_cuda) {
  Vs.setZero();
  Es.setZero();

  if (use_cuda && !aef::is_cuda_initialized()) {
      std::cout << "[HyperfineCalculator] CUDA not initialized but CUDA-based diagonalization called!!" << std::endl;
      std::cout << "[HyperfineCalculator] Converting to Eigen-based diagonalization!!" << std::endl;
      DebugBreak();
      use_cuda = false;
  }

  if (use_cuda) {
      aef::cuda_resize(nBasisElts);
      if (enableDev) {
        aef::diagonalize(H_tot, Es, Vs);
      } else {
        // Free-space: m_f is a good quantum number, want to simultaneously
        // diagonalize H_tot and F_z This is done using the method given in
        // https://math.stackexchange.com/a/4388322 and proven in
        // https://math.stackexchange.com/a/3951339.  However, I'm omitting the
        // randomization portion because it shouldn't be necessary (if there's some
        // small mixing of the m_f it doesn't really matter, and in practice they
        // don't mix.
        constexpr dcomplex t = 100.0; // +15i;
        // naughty hack: H_dev doesn't actually contain anything when enableDev == false, so we can use it
        // as our temporary here instead of making a new temporary matrix
        H_dev = H_tot + t * F_z;
        aef::diagonalize(H_dev, Es, Vs);
        aef::cuda_expectation_values(Vs, H_tot, H_dev);
        Es = H_dev.diagonal();
        H_dev.setZero();
      }
      diagonalized = true;
      return true;
  } else {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
      if (enableDev) {
          // devonshire potential: m_f not a good quantum number --> directly
          // diagonalize H_tot
          solver.compute(H_tot);
          Es = solver.eigenvalues();
          Vs = solver.eigenvectors();
          diagonalized = true;
      } else {
          // Free-space: m_f is a good quantum number, want to simultaneously
          // diagonalize H_tot and F_z This is done using the method given in
          // https://math.stackexchange.com/a/4388322 and proven in
          // https://math.stackexchange.com/a/3951339.  However, I'm omitting the
          // randomization portion because it shouldn't be necessary (if there's some
          // small mixing of the m_f it doesn't really matter, and in practice they
          // don't mix.
          constexpr dcomplex t = 100.0; // +15i;
          Eigen::MatrixXcd temp = H_tot + t * F_z;
          solver.compute(temp);
          Vs = solver.eigenvectors();
          Vst = Vs.adjoint();
#if 0
#define TEST_DIAG
#ifdef TEST_DIAG
          temp = Vst * H_tot * Vs;
          temp.diagonal().setZero();
          assert(temp.isZero(1E-6));

          temp = Vst * H_tot * Vs;
          temp.diagonal().setZero();
          assert(temp.isZero(1E-6));
#endif
#endif
          Es = (Vst * H_tot * Vs).diagonal();
          diagonalized = true;
      }
      return diagonalized;
  }
}

bool HyperfineCalculator::load_matrix_elts(std::string inpath) {
  std::ifstream in(inpath, std::ios::binary);
  return load_matrix_elts(in);
}

bool HyperfineCalculator::save_matrix_elts(std::string outpath) {
  std::ofstream out(outpath,
                    std::ios::binary | std::ios::trunc | std::ios::out);
  return save_matrix_elts(out);
}

bool HyperfineCalculator::load_matrix_elts(fs::path inpath) {
  std::ifstream in(inpath, std::ios::binary);
  return load_matrix_elts(in);
}

bool HyperfineCalculator::save_matrix_elts(fs::path outpath) {
  std::ofstream out(outpath,
                    std::ios::binary | std::ios::trunc | std::ios::out);
  return save_matrix_elts(out);
}

#include <aef/matrix_io.hpp>
#include <algorithm>
#include <iterator>
#include <zstr.hpp>

static constexpr uint8_t MAGIC[8] = {'A', 'e', 'F', 0, 'H', 'D', 'a', 't'};

/*
 * <summary>
 * This enum class describes what each bit of flags means.
 * </summary>
 */
enum hyp_flags : uint16_t {
  // if this flag is set, then the file contains the computed eigenenergies/eigenstates of the Hamiltonian
  FLAG_DIAG = 1,
  // if this flag is set, then the file contains the MDA spherical tensor operators
  FLAG_OKQ = 2,
  // This flag will always be set for files with version >= aefdat_version::rawmat_okq2 so that flags is
  // never all-zero to allow easy distinguishing of files affected by the 4-byte version bug.  Technically
  // this is only neccessary for files with version == aefdat_version::rawmat_okq2 but it may as 
  // note: FLAG_AVOID_ZERO_FLAGS_BUG used to be FLAG_E_DEV_PARMS, however, it was set without the writing
  // code being implemented (oops), so the value of FLAG_E_DEV_PARMS value
  FLAG_AVOID_ZERO_FLAGS_BUG = 4,
  // if this flag is set, then the file contains the specific values of electric field strength and 
  // devonshire coupling constant used to calculate these parts of the Hamiltonian.
  FLAG_E_DEV_PARMS = 8
};

// #define _NO_COMPRESS
uint64_t stream_pos;
bool HyperfineCalculator::load_matrix_elts(std::istream &in) {
#ifndef _NO_COMPRESS
  zstr::istream zin(in);
#else
  auto &zin = in;
#endif
  uint8_t test[8] = {};
  stream_pos = 0;
  std::cout << "READING MAGIC :";
  for (int idx = 0; idx < sizeof(test); idx++) {
    zin >> test[idx];
    stream_pos++;
    std::cout << test[idx];
  }
  std::cout << std::endl;
  // zin.read((char*)test, 8);
  if (memcmp(test, MAGIC, sizeof(MAGIC)) != 0) {
    union {
      uint8_t *u8;
      uint64_t *u64;
    } data;
    data.u8 = test;
    uint64_t rmag = *data.u64;
    data.u8 = (uint8_t *)MAGIC;
    std::cout << "BAD MAGIC " << rmag << " GOOD WOULD BE " << *data.u64
              << std::endl;
    return false;
  } else {
    union {
      uint8_t *u8;
      uint64_t *u64;
    } data;
    data.u8 = test;
    std::cout << "READ CORRECT MAGIC " << *data.u64 << std::endl;
  }
  uint16_t raw_version = 0;
  zin.read((char *)&raw_version, sizeof(raw_version));
  stream_pos += sizeof(raw_version);
  aefdat_version version = static_cast<aefdat_version>(raw_version);

  std::cout << "READ VERSION" <<raw_version << std::endl;

  if (version > CURRENT_VERSION || version < MINIMUM_VERSION) {
    const char *reason =
        (version > CURRENT_VERSION) ? " is too recent." : " is too old.";
    std::cout << "Unsupported version " << static_cast<uint16_t>(version) << reason << std::endl;
    uint16_t flags = 0;
    zin.read((char *)&flags, sizeof(flags));
    std::cout << "FLAGS WOULD BE " << flags << std::endl;
    return false;
  }
  this->load_version = version;
  uint16_t flags;
  zin.read((char *)&flags, sizeof(flags));
  stream_pos += sizeof(flags);

  if (version == aefdat_version::rawmat_okq2 && flags == 0) {
    // hacky bugfix for the "4-byte version bug": the version was sometimes accidentally written as four
    // bytes for part of version aefdat_version::rawmat_okq2. This is because the version was written without
    // being being cast to a u16 first, so sizeof(aefdat_version) bytes woud be written.  GCC and clang chose
    // to use 4 bytes for aefdat_version, so two extra zero bytes were present where the flags should be. 
    // Fortunately,the flags were never written as all-zero during the time this bug was present (and also
    // never will be in newer versions because of FLAG_AVOID_ZERO_FLAGS_BUG). Thus, if "flags" is read as
    // zero in a rawmat_okq2 file, then it was actually the accidental extra two bytes of version, so we need
    // to reread flags to get the correct value -- i.e.:
    // correct: | version - 2B | flags  - 2B | ....
    // old bug: | version - 2B | 0x0000 - 2B | flags - 2B | ....
    zin.read((char *)&flags, sizeof(flags));
  }

  // flags should never be zero for files with version >= aefdat_version::rawmat_okq2
  assert (version <= aefdat_version::rawmat_okq2 || flags != 0);

  std::cout << std::bitset<16>(flags) << " = " << flags << std::endl;

  uint32_t nmax_ = 0;
  zin.read((char *)&nmax_, sizeof(nmax_));
  stream_pos += sizeof(nmax_);
  std::cout << "NMax is " << nmax_ << std::endl;
  nmax = nmax_;
  set_nmax(nmax_);

  // read electric field and devonshire coupling constant
  if (flags & FLAG_E_DEV_PARMS) {
      zin.read((char*)&(this->E_z), sizeof(this->E_z));
      zin.read((char*)&(this->K), sizeof(this->K));
      this->enableDev = (K != 0);
  }

  diagonalized = (flags & FLAG_DIAG);
  // read the various parts of the hamiltonian
  Eigen::read_binary(zin, H_rot.diagonal());
  Eigen::read_binary(zin, H_hfs);
  Eigen::read_binary(zin, H_stk);
  Eigen::read_binary(zin, H_dev);
  Eigen::read_binary(zin, H_tot);

  init = true;

  // if diagonalized, read the eigenenergies and eigenstates
  if (diagonalized) {
    std::cout << "Hamiltonian has been diagonalized" << std::endl;
    Eigen::read_binary(zin, Es);
    Eigen::read_binary(zin, Vs);
  } else {
    std::cout << "Hamiltonian hasn't been diagonalized" << std::endl;
  }

  // dipole-moment spherical tensor operators
  if (version <= aefdat_version::xiff && (flags & FLAG_OKQ)) {
    dkq_init = true;
    Eigen::read_binary(zin, d10);
    Eigen::read_binary(zin, d11);
    Eigen::read_binary(zin, d1t);
  }

  return true;
}

bool HyperfineCalculator::save_matrix_elts(std::ostream &out) {
#ifndef _NO_COMPRESS
  zstr::ostream zout(out);
#else
  auto &zout = out;
#endif
  std::copy(MAGIC, MAGIC + 8, std::ostream_iterator<uint8_t>(zout));
  constexpr aefdat_version version = CURRENT_VERSION;
  uint16_t raw_version = static_cast<uint16_t>(version);
  zout.write((char *)&raw_version, sizeof(raw_version));
  uint16_t flags = 0;

  if (diagonalized) {
    flags |= FLAG_DIAG;
  }

  if (dkq_init) {
    flags |= FLAG_OKQ;
  }
  // always write out the electric field strength and devonshire coupling constant to make LowStateDumper's
  // job easier
  flags |= FLAG_E_DEV_PARMS;
  // flags must never be written out as zero because of a version 3 bug
  flags |= FLAG_AVOID_ZERO_FLAGS_BUG;
  zout.write((char *)&flags, sizeof(flags));
  uint32_t nmax_ = (uint32_t)nmax;
  zout.write((char *)&nmax_, sizeof(nmax_));

  // write out E_z and K
  if (flags & FLAG_E_DEV_PARMS) {
    assert(version >= aefdat_version::rawmat_okq_params);
    zout.write((char*)&(this->E_z), sizeof(this->E_z));
    zout.write((char*)&(this->K), sizeof(this->K));
  }

  // write out the Hamiltonian and its constituent parts
  Eigen::write_binary(zout, H_rot.diagonal());
  Eigen::write_binary(zout, H_hfs);
  Eigen::write_binary(zout, H_stk);
  Eigen::write_binary(zout, H_dev);
  Eigen::write_binary(zout, H_tot);

  // write out Es and Vs if they're initialized
  if (diagonalized) {
    Eigen::write_binary(zout, Es);
    Eigen::write_binary(zout, Vs);
  }

  // write out the spherical tensor operators if they're initialized
  if (flags & FLAG_OKQ) {
    Eigen::write_binary(zout, d10);
    Eigen::write_binary(zout, d11);
    Eigen::write_binary(zout, d1t);
  }

  return true;
}
/*
    aef/HyperfineCalculator.hpp -- contains the HyperfineCalculator class that
    implements the various matrix-element calculations for the 138BaF system.

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
#pragma once
#include "aef.h"
#include <istream>
#include <vector>
#include <filesystem>
#include "Eigen/Eigen"

enum class aefdat_version {
    // original file format: was very similar to format rawmat, but save cod
    invalid = 0,
    // 
    rawmat = 1,
    //
    rawmat_okq = 2,
    // saves more information
    rawmat_okq2 = 3,
    xiff = 4,
    max = xiff
};

class HyperfineCalculator {
public:
    spin nmax;
    size_t nBasisElts;
    bool enableDev;
    aefdat_version load_version;

    // current file format: rawmat_okq2
    static constexpr aefdat_version CURRENT_VERSION = aefdat_version::rawmat_okq2;
    // minimum readable file format: needed because version 0 didn't open files as binary
    static constexpr aefdat_version MINIMUM_VERSION = aefdat_version::rawmat;
private:
    double E_z;
    double K;
    bool init;
    bool diagonalized;
    bool dkq_init;

public:
    std::vector<j_basis_vec> basis;

    Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic> H_rot;
    Eigen::MatrixXcd H_hfs;
    Eigen::MatrixXcd H_stk;
    Eigen::MatrixXcd H_dev;
    Eigen::MatrixXcd F_z;

    Eigen::MatrixXcd H_tot;

    Eigen::VectorXcd Es;
    Eigen::MatrixXcd Vs;
    Eigen::MatrixXcd Vst;

    // dipole-moment matrix operators -- dz = Hstark / muE
    Eigen::MatrixXcd d10;
    Eigen::MatrixXcd d11;
    Eigen::MatrixXcd d1t;

public:
    /// <summary>
    /// Constructor
    /// </summary>
    /// <param name="nmax_">Maximum value of the n quantum number.  There will be 4*nmax*nmax states in the basis</param>
    /// <param name="E_z">External electric field strength to use for Stark calculations (MHz/D).</param>
    /// <param name="K">Devonshire coupling constant (MHz)</param>
    HyperfineCalculator(spin nmax_ = 0.0, double E_z = 1.0, double K=0.0);
    ~HyperfineCalculator();

    bool isDiagonalized() {
        return diagonalized;
    }
    bool calculate_matrix_elts();
    bool diagonalize_H();

    /// <summary>
    /// Are the orientation spherical tensor operatorsinitialized
    /// </summary>
    /// <returns></returns>
    bool isDkqInit() { return dkq_init;}

    /// <summary>
    /// Calculate the dipole-m
    /// </summary>
    /// <returns></returns>
    bool calc_dkq();

    bool load_matrix_elts(std::string inpath);
    bool save_matrix_elts(std::string out);

    bool load_matrix_elts(std::filesystem::path inpath);
    bool save_matrix_elts(std::filesystem::path out);

    bool load_matrix_elts(std::istream& in);
    /// <summary>
    /// Save the various operators to an AeF0Dat stream.  Only supports writing as the latest version
    /// </summary>
    /// <param name="out"></param>
    /// <returns></returns>
    bool save_matrix_elts(std::ostream& out);

    void set_nmax(spin nmax_);

    inline dcomplex eval_H(Eigen::VectorXcd& v, Eigen::VectorXcd& w) {
        return (v.transpose() * H_tot * w)(0, 0);
    }
};

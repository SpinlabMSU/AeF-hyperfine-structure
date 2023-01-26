#pragma once
#include "aef.h"
#include <istream>
#include <vector>
#include "Eigen/Eigen"

class HyperfineCalculator {
public:
    spin nmax;
    size_t nBasisElts;
    bool enableDev;
private:
    double E_z;
    bool init;
    bool diagonalized;

public:
    std::vector<j_basis_vec> basis;

    Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic> H_rot;
    Eigen::MatrixXcd H_hfs;
    Eigen::MatrixXcd H_stk;
    Eigen::MatrixXcd H_dev;

    Eigen::MatrixXcd H_tot;

    Eigen::VectorXcd Es;
    Eigen::MatrixXcd Vs;

public:
    HyperfineCalculator(spin nmax_ = 0.0, double E_z = 1.0, bool enableDev = false);
    ~HyperfineCalculator();

    bool calculate_matrix_elts();
    bool diagonalize_H();
    bool load_matrix_elts(std::istream& in);
    bool save_matrix_elts(std::ostream& out);

    void set_nmax(spin nmax_);

    inline dcomplex eval_H(Eigen::VectorXcd& v, Eigen::VectorXcd& w) {
        return (v.transpose() * H_tot * w)(0, 0);
    }
};

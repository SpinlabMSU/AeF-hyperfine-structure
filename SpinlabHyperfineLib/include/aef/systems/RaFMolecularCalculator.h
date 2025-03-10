#pragma once
#include "aef/MolecularSystem.h"
#include <aef/j_basis_vec.h>
namespace aef {
    class RaFMolecularCalculator : public IMolecularCalculator {
    private:
        spin nmax;
        size_t nBasisElts;
    public:
        std::vector<jf_basis_vec> basis;
        RaFMolecularCalculator();
        RaFMolecularCalculator(spin nmax_);
        virtual ~RaFMolecularCalculator();


        virtual ResultCode get_parameter(std::string id, double &out);
        virtual ResultCode set_parameter(std::string id, double value);
        virtual spin get_nmax();
        virtual void set_nmax(spin nmax);
        virtual size_t get_nBasisElts() {
            return nBasisElts;
        };

        // hamiltonian
        virtual ResultCode calculate_H_rot(Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic>& H);
        virtual ResultCode calculate_H_hfs(Eigen::MatrixXcd& H);
        virtual ResultCode calculate_H_dev(Eigen::MatrixXcd& H);
        virtual ResultCode calculate_H_stk(Eigen::MatrixXcd& H);


        // operators
        virtual void calculate_F_z(Eigen::MatrixXcd& F_z);
        virtual ResultCode calculate_dkq(Eigen::MatrixXcd& d, int q);
        virtual void calculate_d1t(Eigen::MatrixXcd& H);
        virtual void calculate_d10(Eigen::MatrixXcd& H);
        virtual void calculate_d11(Eigen::MatrixXcd& H);

        // IO
        virtual void load(std::istream& in);
        virtual void save(std::ostream& out);
    };

}
#pragma once
#include "aef/MolecularSystem.h"
#include <aef/j_basis_vec.h>
namespace aef {
    class BaFMolecularCalculator : public IMolecularCalculator {
    private:
        spin nmax;
        size_t nBasisElts;
        std::vector<universal_diatomic_basis_vec> lowest_states;
    public:
        std::vector<j_basis_vec> basis;
        BaFMolecularCalculator();
        BaFMolecularCalculator(spin nmax_);
        virtual ~BaFMolecularCalculator();


        virtual ResultCode get_parameter(std::string id, double &out);
        virtual ResultCode set_parameter(std::string id, double value);
        virtual spin get_nmax();
        virtual void set_nmax(spin nmax);
        virtual universal_diatomic_basis_vec get_basis_ket(int idx);
        virtual int get_index(universal_diatomic_basis_vec v);

        virtual size_t get_nBasisElts() {
            return nBasisElts;
        };

        // hamiltonian
        virtual ResultCode calculate_H_rot(Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic>& H);
        virtual ResultCode calculate_H_hfs(Eigen::MatrixXcd& H);
        virtual ResultCode calculate_H_dev(Eigen::MatrixXcd& H, double K=1.0);
        virtual ResultCode calculate_H_stk(Eigen::MatrixXcd& H, double E_z=1.0);


        // operators
        virtual void calculate_F_z(Eigen::MatrixXcd& F_z);
        virtual ResultCode calculate_dkq(Eigen::MatrixXcd& d, int q);
        virtual void calculate_d1t(Eigen::MatrixXcd& H);
        virtual void calculate_d10(Eigen::MatrixXcd& H);
        virtual void calculate_d11(Eigen::MatrixXcd& H);

        virtual std::array<dcomplex, 3> molec_edm(int kdx1, int kdx2);
        virtual std::array<dcomplex, 3> molec_mdm(int kdx1, int kdx2);
        virtual void calculate_S_dot_ina(Eigen::MatrixXcd& A);
        virtual void calculate_I1_dot_ina(Eigen::MatrixXcd& A);
        virtual void calculate_I2_dot_ina(Eigen::MatrixXcd& A);

        // IO
        virtual void load(std::istream& in);
        virtual void save(std::ostream& out);
        virtual const char* get_calc_type();

        virtual int get_num_orientations() {
            return 6;
        }
        virtual int get_lowest_rotational_state_size();

        virtual std::vector<universal_diatomic_basis_vec> get_lowest_states();
    };

}
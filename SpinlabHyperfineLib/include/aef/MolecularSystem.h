/*
    aef/MolecularSystem.h -- class that implements matrix-element calculations
    for the 225RaF molecule.

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
#ifndef _AEF_MOLECULAR_SYSTEM_H

#define _AEF_MOLECULAR_SYSTEM_H 1
#include "Eigen/Eigen"
#include "aef.h"
#include "jf_basis_vec.h"
#include <filesystem>
#include <istream>
#include <unordered_map>
#include <span>
#include <unordered_set>
#include <vector>

#include <aef/operators/IOperator.h>


namespace aef {

    /// <summary>
    /// This struct has information
    /// 
    /// Under normal circumstances, 
    /// </summary>
    struct universal_diatomic_basis_vec {
        /// <summary>
        /// The quantum number describing the magnitude of nuclues #1's spin (the nucleus with the larger hyperfine structure).
        /// 
        /// AeF-hyperfine-structure supports $I_1$ = 1/2.
        /// </summary>
        spin i1;
        /// <summary>
        /// The quantum number describing the magnitude of nucleus #2's spin (the nucleus with the smaller hyperfine structure).
        /// 
        /// AeF-hyperfine-structure supports $I_2$ = 0, 1/2.
        /// </summary>
        spin i2;
        /// <summary>
        /// The quantum number describing the magnitude of the molecule's net electron spin.
        /// 
        /// For electronic states AeF-hyperfine-structure currently supports (^2\Sigma^+) this is 1/2.
        /// </summary>
        spin s;
        /// <summary>
        /// The quantum number describing the magnitude of the the molecule's net electron orbital angular momentum.
        /// 
        /// For electronic states AeF-hyperfine-structure currently supports (^2\Sigma^+) this is 0.
        /// </summary>
        spin l;
        /// <summary>
        /// The quantum number describing the magnitude of the "R" internuclear-orbital angular momentum.
        /// 
        /// Since AeF-hyperfine-structure
        /// </summary>
        spin r;
        /// <summary>
        /// The quantum number describing the magnitude of the "N" rotational angular momentum (technically R
        /// </summary>
        spin n;
        /// <summary>
        /// The quantum number describing the magnitude of the "J" total non-nuclear-spin angular momentum
        /// </summary>
        spin j;
        /// <summary>
        /// The quantum number describing the magnitude of the "F_1" angular momentum including one nuclear spin
        /// </summary>
        spin f_1;
        /// <summary>
        /// The quantum number describing the magnitude of "F", the total angular momentum
        /// </summary>
        spin f;
        /// <summary>
        /// The "magnetic" quantum number for the total angular mometum.
        /// </summary>
        spin m_f;

        enum class coupling_type {
            invalid,
            j_basis,
            jf_basis
        };

        coupling_type type;

        universal_diatomic_basis_vec(): i1(0), i2(0), l(0), s(0), r(0), f_1(0), n(0), j(0), f(0), m_f(0) {
            type = coupling_type::invalid;
        }
        universal_diatomic_basis_vec(spin n_, spin j_, spin f_, spin m_f_): i1(half), i2(0), l(0), s(half), r(n_), f_1(f_),
            n(n_), j(j_), f(f_), m_f(m_f_){
            type = coupling_type::j_basis;
        }
        universal_diatomic_basis_vec(spin n_, spin j_, spin f_1_, spin f_, spin m_f_) :i1(half), i2(0), l(0), s(half), r(n_),
            n(n_), j(j_), f_1(f_1_), f(f_), m_f(m_f_) {
            type = coupling_type::jf_basis;
        }

        /// <summary>
        /// Descibes this state as a ket
        /// </summary>
        /// <returns>A string</returns>
        std::string ket_string();
        /// <summary>
        /// Returns a string that can be used in a CSV field without quoting
        /// </summary>
        /// <returns>description suitable for CSV</returns>
        std::string ket_csv_str(); 
    };




    constexpr size_t max_op_id_len = 256;
    class IMolecularCalculator {
    public:
        virtual ResultCode get_parameter(std::string id, double& out) = 0;
        virtual ResultCode set_parameter(std::string id, double value) = 0;
        virtual void set_nmax(spin nmax) = 0;
        virtual size_t get_nBasisElts() = 0;
        virtual universal_diatomic_basis_vec get_basis_ket(int idx) = 0;
        virtual int get_index(universal_diatomic_basis_vec v) = 0;

        // hamiltonian
        virtual ResultCode calculate_H_rot(Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic>& H) = 0;
        virtual ResultCode calculate_H_hfs(Eigen::MatrixXcd& H) = 0;
        virtual ResultCode calculate_H_dev(Eigen::MatrixXcd& H, double K=1) = 0;
        virtual ResultCode calculate_H_stk(Eigen::MatrixXcd& H, double E_z=1) = 0;
        
        
        // operators
        virtual void calculate_F_z(Eigen::MatrixXcd& F_z) = 0;
        virtual ResultCode calculate_dkq(Eigen::MatrixXcd& d, int q) = 0;
        virtual void calculate_d1t(Eigen::MatrixXcd& A) = 0;
        virtual void calculate_d10(Eigen::MatrixXcd& A) = 0;
        virtual void calculate_d11(Eigen::MatrixXcd& A) = 0;
        virtual void calculate_S_dot_ina(Eigen::MatrixXcd&A) = 0;
        virtual void calculate_I1_dot_ina(Eigen::MatrixXcd& A) = 0;
        virtual void calculate_I2_dot_ina(Eigen::MatrixXcd& A) = 0;

        // load/save any state
        virtual void load(std::istream& in) = 0;
        virtual void save(std::ostream& out) = 0;
        virtual const char* get_calc_type() = 0;

        virtual int get_num_orientations() = 0;
        virtual int get_lowest_rotational_state_size() = 0;
        virtual std::vector<universal_diatomic_basis_vec> get_lowest_states() = 0;
        
        virtual std::array<dcomplex, 3> molec_edm(int kdx1, int kdx2) = 0;
        virtual std::array<dcomplex, 3> molec_mdm(int kdx1, int kdx2) = 0;
    public:
        // convenience function for use as a molcalcmaker
        template<class T> static IMolecularCalculator* createInstance() {
            return new T;
        }
        typedef IMolecularCalculator* (*pfnMolCalcMaker)();
        static void registerMolCalcType(std::string name, pfnMolCalcMaker ctor);
        static IMolecularCalculator* makeCalculatorOfType(std::string name);
        static void register_default_types();
    };

#include "io/molsys_io.h"

    /// <summary>
    /// The aef::MolecularSystem class implements a generic
    /// </summary>
    class MolecularSystem {
    private:
        std::unordered_map<std::string, Eigen::MatrixXcd*> opMatMap;
        std::unordered_map<std::string, aef::operators::IOperator*> opMap;
    public:
        spin nmax;
        size_t nBasisElts;
    private:
        IMolecularCalculator *calc;
        bool init;
        bool diagonalized;
        bool dkq_init;
    public:
        double E_z;
        double K;
        bool enableDev;

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
        /// <summary>
        /// Use this constructor to 
        /// </summary>
        /// <param name="calc"></param>
        /// <param name="nmax_"></param>
        /// <param name="E_z_"></param>
        /// <param name="K"></param>
        MolecularSystem(IMolecularCalculator* calc, spin nmax_, double E_z_ = 0, double K=0.0);
        /// <summary>
        /// Use this constructor to make molecularsystems that you will immediately call load on.
        /// No operations other than load can be performed on a MolecularSystem constructed by this
        /// constructor until a load succeeds.
        /// </summary>
        MolecularSystem();
        ~MolecularSystem();

        IMolecularCalculator* get_calc();

        void set_nmax(spin nmax_); // 
        void calculate_matrix_elts();
        aef::ResultCode diagonalize(); // always uses aef::matrix now


        aef::ResultCode load(std::istream& in, const char *path=nullptr);
        aef::ResultCode save(std::ostream& out, const char *path=nullptr);

        // convienence methods
        aef::ResultCode load(std::string inpath);
        aef::ResultCode save(std::string out);

        aef::ResultCode load(std::filesystem::path inpath);
        aef::ResultCode save(std::filesystem::path out);

        inline dcomplex eval_H(Eigen::VectorXcd& v, Eigen::VectorXcd& w) {
            return (v.transpose() * H_tot * w)(0, 0);
        }

    private:
        aef::ResultCode write_chunk(std::ostream& out, void* chdr, void* data);
        aef::ResultCode read_chunk(std::istream &in, void *chdr, void *dst);

        aef::ResultCode write_matrix(std::ostream& out, Eigen::MatrixXcd* mat, uint32_t matnam);
        aef::ResultCode write_vector(std::ostream& out, Eigen::VectorXcd* vec, uint32_t matnam);

    /// <summary>
    /// This section contains code for PTFW2
    /// </summary>
    public:
        /// <summary>
        /// Finds operator with given ID
        /// </summary>
        /// <param name="id"></param>
        /// <returns></returns>
        aef::operators::IOperator* getOperator(const std::string& id);
        bool hasOperator(const std::string& id) {
            return nullptr != getOperator(id);
        }

        /// <summary>
        /// Registers an operator with the given id
        /// </summary>
        /// <param name="id"></param>
        /// <param name="op"></param>
        void addOperator(const std::string& id, aef::operators::IOperator* op);

        /// <summary>
        /// Evaluate the matrix forms of all registered operators
        /// </summary>
        void evaluate(void);

        /// <summary>
        /// Get the matrix form of the given operator in the default basis
        /// </summary>
        /// <param name="id"></param>
        /// <returns></returns>
        Eigen::MatrixXcd* getOperatorMatrix(const std::string& id);

        /// <summary>
        /// Gets the matrix element of the given operator in the energy-eigenstate basis
        /// </summary>
        /// <param name="eidx1"></param>
        /// <param name="eidx2"></param>
        /// <returns></returns>
        dcomplex get_matrix_element(const std::string& id, int eidx1, int eidx2);

        /// <summary>
        /// Return the expectation value of the given operator using the
        /// </summary>
        /// <param name="id"></param>
        /// <param name="eidx1"></param>
        /// <returns></returns>
        dcomplex expectation_value(const std::string& id, int eidx1);

        /// <summary>
        /// Perturbatively calculate the leading order changes in energy of the energy eigenstates resulting
        /// from the given operator.
        /// 
        /// Limitation: currently only handles first order.
        /// </summary>
        /// <param name="id">operator ID</param>
        /// <param name="output"></param>
        /// <param name="workspace"></param>
        aef::ResultCode delta_E_lo(const std::string& id, Eigen::VectorXcd& output, Eigen::MatrixXcd* workspace = nullptr);
    };
};

using aef::universal_diatomic_basis_vec;
std::ostream& operator<<(std::ostream& os, universal_diatomic_basis_vec& v);
bool operator == (const universal_diatomic_basis_vec& v1, const universal_diatomic_basis_vec& v2);

template <> struct fmt::formatter<universal_diatomic_basis_vec> : fmt::formatter<std::string> {
    auto format(universal_diatomic_basis_vec v, format_context& ctx) const {
        if (v.type == universal_diatomic_basis_vec::coupling_type::j_basis) {
            return formatter<std::string>::format(fmt::format("|n={},j={},f={},m_f={}>", v.n, v.j, v.f, v.m_f), ctx);
        } else if (v.type == universal_diatomic_basis_vec::coupling_type::jf_basis) {
            return formatter<std::string>::format(fmt::format("|n={},j={},f_1={},f={},m_f={}>", v.n, v.j, v.f_1, v.f, v.m_f), ctx);
        } else {
            return formatter<std::string>::format(fmt::format("|i_1={}, i_2={}, s={}, l={}, r={}, n={},j={},f_1={},f={},m_f={}>", 
                v.i1, v.i2, v.s, v.l, v.r,v.n, v.j, v.f_1, v.f, v.m_f), ctx);
        }
    }
};
#endif //_AEF_MOLECULAR_SYSTEM_H
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

    class IMolecularCalculator {
    public:
        virtual ResultCode get_parameter(std::string id, double& out) = 0;
        virtual ResultCode set_parameter(std::string id, double value) = 0;
        virtual void set_nmax(spin nmax) = 0;
        virtual size_t get_nBasisElts() = 0;

        // hamiltonian
        virtual ResultCode calculate_H_rot(Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic>& H) = 0;
        virtual ResultCode calculate_H_hfs(Eigen::MatrixXcd& H) = 0;
        virtual ResultCode calculate_H_dev(Eigen::MatrixXcd& H) = 0;
        virtual ResultCode calculate_H_stk(Eigen::MatrixXcd& H) = 0;
        
        
        // operators
        virtual void calculate_F_z(Eigen::MatrixXcd& F_z) = 0;
        virtual ResultCode calculate_dkq(Eigen::MatrixXcd& d, int q) = 0;
        virtual void calculate_d1t(Eigen::MatrixXcd& H) = 0;
        virtual void calculate_d10(Eigen::MatrixXcd& H) = 0;
        virtual void calculate_d11(Eigen::MatrixXcd& H) = 0;
        

    public:
        // convenience function for use as a molcalcmaker
        template<class T> static IMolecularCalculator* createInstance() {
            return new T;
        }
        typedef IMolecularCalculator* (*pfnMolCalcMaker)();
        static void registerMolCalcType(std::string name, pfnMolCalcMaker ctor);
        static IMolecularCalculator* makeCalculatorOfType(std::string name);
    };

    struct ExternalFieldParameters {
        double E_z;
        double K;
    };

    /// <summary>
    /// The aef::MolecularSystem class implements a generic
    /// </summary>
    class MolecularSystem {
    private:
        std::unordered_map<std::string, Eigen::MatrixXcd*> opMatMap;
        std::unordered_map<std::string, aef::operators::IOperator*> opMap;

        spin nmax;
        size_t nBasisElts;
        bool enableDev;
    private:
        IMolecularCalculator &calc;
        bool init;
        bool diagonalized;
        bool dkq_init;
    public:
        double E_z;
        double K;

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

        MolecularSystem(IMolecularCalculator& calc, spin nmax_, double E_z_ = 0, double K=0.0);
        ~MolecularSystem();

        void set_nmax(spin nmax_); // 
        void calculate_matrix_elts();
        void calculate_dkq(); // 
        aef::ResultCode diagonalize(); // always uses aef::matrix now


        aef::ResultCode load(std::istream& in);
        aef::ResultCode save(std::ostream& out);

        // convienence methods
        aef::ResultCode load(std::string inpath);
        aef::ResultCode save(std::string out);

        aef::ResultCode load(std::filesystem::path inpath);
        aef::ResultCode save(std::filesystem::path out);

        inline dcomplex eval_H(Eigen::VectorXcd& v, Eigen::VectorXcd& w) {
            return (v.transpose() * H_tot * w)(0, 0);
        }


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
#endif //_AEF_MOLECULAR_SYSTEM_H
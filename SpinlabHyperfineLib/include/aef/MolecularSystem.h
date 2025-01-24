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
#include <vector>

#include <aef/operators/IOperator.h>


namespace aef {

    class IMolecularCalculator {
    public:
        virtual size_t nBasisElts() = 0;

        // hamiltonian
        virtual ResultCode calculate_H_rot(Eigen::DiagonalMatrix<dcomplex, Eigen::Dynamic>& H) = 0;
        virtual ResultCode calculate_H_hfs(Eigen::MatrixXcd& H) = 0;
        virtual ResultCode calculate_H_dev(Eigen::MatrixXcd& H) = 0;
        virtual ResultCode calculate_H_stk(Eigen::MatrixXcd& H) = 0;
        
        
        // operators
        virtual ResultCode calculate_dkq(Eigen::MatrixXcd& d, int q) = 0;
        virtual void calculate_d10(Eigen::MatrixXcd& H) = 0;

        //
        virtual ResultCode getOperatorList(IOperator<)

    public:
        template<class T> static IMolecularCalculator* createInstance() {
            return new T;
        }
        typedef IMolecularCalculator* (*pfnMolCalcMaker)();
        static void registerMolCalcType(std::string name, pfnMolCalcMaker ctor);
        static IMolecularCalculator *makeCalculatorOfType
    private:


    };

    struct ExternalFieldParameters {
        double E_z;
        double K;
    };

    /// <summary>
    /// The aef::MolecularSystem class implements a generic
    /// </summary>
    class MolecularSystem {
        // class ifnfo
        typedef IMolecularCalculator(*molCalcMaker)();
        static std::unordered_map<std::string, molCalcMaker> molCalcMakerMap;
        // 
    private:
        std::unordered_map<std::string, Eigen::MatrixXcd*> operators;

        spin nmax;
        size_t nBasisElts;
        bool enableDev;
    private:
        IMolecularCalculator *calc;
        bool init;
        bool diagonalized;
        bool dkq_init;
    public:

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

        MolecularSystem(IMolecularCalculator& calc);
        ~MolecularSystem();



        void set_nmax(spin nmax_);
        void calculate_matrix_elts();
        void calculate_dkq();



        bool load(std::istream& in);
        bool save(std::ostream& out);

        // convienence methods
        bool load_matrix_elts(std::string inpath);
        bool save_matrix_elts(std::string out);

        bool load_matrix_elts(std::filesystem::path inpath);
        bool save_matrix_elts(std::filesystem::path out);


        inline dcomplex eval_H(Eigen::VectorXcd& v, Eigen::VectorXcd& w) {
            return (v.transpose() * H_tot * w)(0, 0);
        }
    };
};
#endif //_AEF_MOLECULAR_SYSTEM_H
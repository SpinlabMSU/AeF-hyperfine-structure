/*
    aef/operators/PerturbationFramework.h -- contains various utility functions

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
#ifndef _AEF_OPERATORS_PERTURBATIONFRAMEWORK_H
#define _AEF_OPERATORS_PERTURBATIONFRAMEWORK_H 1

#include <aef/aef.h>
#include <aef/operators/IOperator.h>

namespace aef::operators{
    class PerturbationFramework {
        using basis_vec = aef::j_basis_vec;
        using System = HyperfineCalculator;

        System* sys;
        std::unordered_map<std::string, IKetOperator<basis_vec>*> opMap;
        std::unordered_map<std::string, Eigen::MatrixXcd*> opMatMap;
        void* pImpl;
    // eventual interface
    public:
        PerturbationFramework(System *sys_);
        ~PerturbationFramework();
        System* get_system();
        void set_basis_system(System* calc);

        IKetOperator<basis_vec>* getOperator(const std::string& id);
        bool hasOperator(const std::string& id) {
            return nullptr != getOperator(id);
        }
        void addOperator(const std::string &id, IKetOperator<basis_vec>* op);

        /// <summary>
        /// Perform
        /// </summary>
        void evaluate(void);

        Eigen::MatrixXcd* getOperatorMatrix(const std::string& id);
    
        /// <summary>
        /// Gets the matrix element
        /// </summary>
        /// <param name="eidx1"></param>
        /// <param name="eidx2"></param>
        /// <returns></returns>
        dcomplex get_matrix_element(const std::string &id, int eidx1, int eidx2);
        /// <summary>
        /// Return the expectation value    
        /// </summary>
        /// <param name="id"></param>
        /// <param name="eidx1"></param>
        /// <returns></returns>
        dcomplex expectation_value(const std::string &id, int eidx1);

        /// <summary>
        /// Calculate the leading order change
        /// 
        /// Limitation: currently only handles first order.
        /// </summary>
        /// <param name="id"></param>
        /// <param name="output"></param>
        /// <param name="workspace"></param>
        aef::ResultCode delta_E_lo(const std::string &id, Eigen::VectorXcd &output, Eigen::MatrixXcd *workspace=nullptr);
    };
}
#endif //_AEF_OPERATORS_PERTURBATIONFRAMEWORK_H
/*
    aef/operators/IOperator.h -- contains various utility functions

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
#ifndef _AEF_OPERATOR_IOPERATOR_H
#define _AEF_OPERATOR_IOPERATOR_H 1
#pragma once

#include <aef/aef.h>
#include <aef/operators/basis_ket.h>

namespace aef::operators {

    struct OperatorInfo {
        std::string name;
        std::string description;
        uint32_t is_hermitian : 1;
        uint32_t is_potential : 1;
    };

    /// <summary>
    /// IOperator&lt;BasisKet&gt;
    /// </summary>
    /// <typeparam name="basis_ket">The Basis Ket type</typeparam>
    template <IBasisKet<int> basis_ket> class IOperator {
    public:
        virtual dcomplex matrixElement(basis_ket k1, basis_ket k2) = 0;
        virtual void fillMatrix(Eigen::SparseMatrix<dcomplex> &matrix) = 0;
        virtual void fillMatrix(Eigen::MatrixXcd& matrix) = 0;
    
        virtual OperatorInfo *getInfo() = 0;

        IOperator() {}
        //
        // as a -
        virtual ~IOperator() {}
    };
}
#endif //_AEF_OPERATOR_IOPERATOR_H
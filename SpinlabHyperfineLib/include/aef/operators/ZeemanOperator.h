/*
    aef/operators/ZeemanOperator.h -- contains various utility functions

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

#ifndef _AEF_OPERATOR_ZEEMANOPERATOR_H
#define _AEF_OPERATOR_ZEEMANOPERATOR_H 1
#pragma once

#include <aef/aef.h>
#include "IOperator.h"

namespace aef::operators {
    class ZeemanOperator : public IOperator<aef::j_basis_vec> {
        OperatorInfo info;
        double B;
        union {
            std::array<double, 3> B_vec;
            struct {
                double B_x;
                double B_y;
                double B_z;
            };
        };
    public:
        ZeemanOperator(double B=1.0);
        ZeemanOperator(std::array<double, 3> dir);
        ~ZeemanOperator();

        using basis_ket = aef::j_basis_vec;
        virtual dcomplex matrixElement(basis_ket k1, basis_ket k2);
        virtual void fillMatrix(Eigen::SparseMatrix<dcomplex>& matrix);
        virtual void fillMatrix(Eigen::MatrixXcd& matrix);

        virtual OperatorInfo* getInfo();

    };
};

#endif //_AEF_OPERATOR_ZEEMANOPERATOR_H
#ifndef _AEF_OPERATOR_IOPERATOR_H
#define _AEF_OPERATOR_IOPERATOR_H 1
#pragma once

#include <aef/aef.h>
#include <aef/operators/basis_ket.h>

namespace aef::operators {

    struct OperatorInfo {
        char* name;
        uint32_t is_hermitian : 1;
    };

    /// <summary>
    /// IOperator&lt;BasisKet&gt;
    /// </summary>
    /// <typeparam name="basis_ket">The Basis Ket type</typeparam>
    template <IBasisKet<int> basis_ket> class IOperator {
        virtual dcomplex matrixElement(basis_ket k1, basis_ket k2) = 0;
        virtual void fillMatrix(Eigen::SparseMatrix<dcomplex> &matrix) = 0;
        virtual void fillMatrix(Eigen::MatrixXcd& matrix) = 0;
    
        virtual OperatorInfo getInfo() = 0;

        //
        // as a -
        virtual ~IOperator() {}
    };
}
#endif //_AEF_OPERATOR_IOPERATOR_H
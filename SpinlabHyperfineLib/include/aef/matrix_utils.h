/*
    aef/matrix_utils.h -- provides some extra utility functions for matricies.

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
#ifndef _AEF_MATRIX_UTILS
#define _AEF_MATRIX_UTILS

#include "aef_types.h"
#ifndef __NVCC__
#include "aef_utils.h"
#endif
#include "gcem.hpp"
#include <complex>
#include <numbers>
#include <Eigen/Eigen>
#include <aef/MatrixOpBackend.h>

namespace aef::matrix {
    ResultCode init(BackendType hint, int argc, char** argv);
    IMatrixOpBackend *get_backend();
    // always
    IMatrixOpBackend *get_fallback_backend();

    /// <summary>
    /// Close out the backend
    /// </summary>
    /// <returns>Result code</returns>
    ResultCode shutdown();

    /// <summary>
    /// Allocate enough space for 
    /// </summary>
    /// <param name="nMaxDim">Maximum matrix dimension</param>
    /// <returns>Result code</returns>
    inline ResultCode set_max_size(int nMaxDim) {
        return get_backend()->set_max_size(nMaxDim);
    }

    /// <summary>
    /// Multiplies two matricies
    /// </summary>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <param name="out"></param>
    /// <returns>Result code</returns>
    inline ResultCode multiply(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
        return get_backend()->multiply(A, B, out);
    }

    /// <summary>
    /// Computes the commutator [A, B]
    /// </summary>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <param name="out"></param>
    /// <returns></returns>
    inline ResultCode commutator(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
        return get_backend()->commutator(A, B, out);
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="out"></param>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <param name="workspace"></param>
    /// <returns></returns>
    ResultCode commutes(bool &out, Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd *workspace=nullptr, double prec=1e-8);
    bool commutes(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd *workspace=nullptr, double prec = 1e-8);

    /// <summary>
    /// Computes the "group action" UAU^{-1} of 
    /// </summary>
    /// <param name="U">A unitary matrix</param>
    /// <param name="A">A general matrix</param>
    /// <param name="out">The output matrix</param>
    /// <returns>Result code</returns>
    inline ResultCode group_action(Eigen::MatrixXcd& out, Eigen::MatrixXcd& U, Eigen::MatrixXcd& A) {
        return get_backend()->group_action(out, U, A);
    }
    /// <summary>
    /// Calculates the expectation value of operator A in state v1
    /// </summary>
    /// <param name="out">expectation value</param>
    /// <param name="v1"></param>
    /// <param name="A"></param>
    /// <returns></returns>
    inline ResultCode expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A) {
        return get_backend()->expectation_value(out, v1, A);
    }
    /// <summary>
    /// Calculates the matrix element <v
    /// </summary>
    /// <param name="out"></param>
    /// <param name="v1"></param>
    /// <param name="A"></param>
    /// <param name="v2"></param>
    /// <returns></returns>
    inline ResultCode matrix_element(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::VectorXcd& v2) {
        return get_backend()->matrix_element(out, v1, A, v2);
    }

    /// <summary>
    /// Diagonalizes complex Hermitian matricies (ZHEEV)
    /// </summary>
    /// <param name="mat">The matrix to diagonalize, must be hermitian, not overwiten</param>
    /// <param name="evals"></param>
    /// <param name="evecs"></param>
    inline ResultCode diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) {
        return get_backend()->diagonalize(mat, evals, evecs);
    }
};


namespace aef {
    /// <summary>
    /// Simultaneously diagonalizes matricies A and B under the following two
    /// assumptions:
    ///     1) A and B are diagonalizable.
    ///     2) [A, B] = 0 (A and B must commute
    ///
    /// Uses the algorithm provided by https://math.stackexchange.com/a/4388322:
    /// with probability 1, choosing random $t\in\mathbb{R}$ and diagonalizing A+tB
    /// will provide an eigenbasis of both
    ///
    ///
    /// </summary>
    /// <typeparam name="Matrix"></typeparam>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <param name="U"></param>
    /// <param name="Ut"></param>
    template <class Matrix>
    void simultaneously_diagonalize(Matrix& U, Matrix& Ut, const Matrix& A,
        const Matrix& B, const Matrix& C) {
        //
#ifdef MATRIX_DEBUG
  // ensure
#endif
        typedef typename Matrix::Scalar Scalar;
        Scalar t = 0.0;
        Eigen::SelfAdjointEigenSolver<Matrix> solver;
        constexpr double prec = 1e-4;
        while (true) {
            // randomly choose t --> 0 not in range, also,
            t = genrandom(0.5, 100.0);
            Matrix AtB = A + t * B;
            solver.compute(AtB);

            U = solver.eigenvectors();
            Ut = U.adjoint();

            // check if A diagonalized
            Matrix check = Ut * A * U;
            check.diagonal().setZero();
            bool agood = check.isZero(prec);

            // check if B diagonalized
            check = Ut * B * U;
            check.diagonal().setZero();
            bool bgood = check.isZero(prec);

            // if both are diagonalized then return
            if (agood && bgood)
                return;
            //std::cout << fmt::format("{} FAIL FAIL FAIL", t) << std::endl;
            //DebugBreak();
        }
    }

    template <class Matrix> bool commutes(Matrix& A, Matrix& B) {
        auto comm = A * B - B * A;
        return comm.isZero(1E-8);
    }

#ifdef __INTELLISENSE__
    static void foo() {
        Eigen::MatrixXcd d;
        simultaneously_diagonalize(d, d, d, d, d);
        bool b = commutes(d, d);
    }
#endif
}
#endif
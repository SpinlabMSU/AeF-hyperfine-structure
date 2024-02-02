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

namespace aef {
    /// <summary>
    /// Initializes cuda
    /// </summary>
    /// <param name="argc">Argument count</param>
    /// <param name="argv">Argument list pointer</param>
    /// <returns>True if initializing cuda succeded, false otherwise (note: as currently implemented, actually just crashes)</returns>
    bool init_cuda(int argc, const char** argv);
    /// <summary>
    /// Shuts down cuda
    /// </summary>
    /// <returns></returns>
    bool shutdown_cuda();
    
    bool is_cuda_initialized();

    /// <summary>
    /// Sets the size of the GPU-side CUDA buffers
    /// </summary>
    /// <param name="n"></param>
    void cuda_resize(int n);

    void log_dev_props_info(std::ostream& out);

    dcomplex cuda_expectation_value(Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::MatrixXcd& v2);
    dcomplex cuda_expectation_value(Eigen::VectorXcd& v1, Eigen::MatrixXcd& A);
    void cuda_expectation_values(Eigen::MatrixXcd& U, Eigen::MatrixXcd& A, Eigen::MatrixXcd &out);

    /// <summary>
    /// 
    /// </summary>
    /// <param name="mat"></param>
    /// <param name="evals"></param>
    /// <param name="evecs"></param>
    void diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs);

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
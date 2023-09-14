#pragma once
#ifndef _AEF_MATRIX_UTILS
#define _AEF_MATRIX_UTILS

#include "aef_types.h"
#include "gcem.hpp"
#include "gsl/gsl_sf_coupling.h"
#include <format>
#include <complex>
#include <numbers>
#include <random>
#include <Eigen/Eigen>
#include "pcg/pcg_random.hpp"
#include <type_traits>

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
void simultaneously_diagonalize(Matrix &U, Matrix &Ut, const Matrix &A,
                                const Matrix &B, const Matrix &C) {
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
    std::cout << std::format("{} FAIL FAIL FAIL", t) << std::endl;
    DebugBreak();
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

#endif
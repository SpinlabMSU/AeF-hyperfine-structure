/*
    aef/aef_utils.h -- contains various utility functions
    
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
#ifndef _AEF_UTILS_H
#define _AEF_UTILS_H 1
#pragma once

#include "aef_types.h"
#include "gcem.hpp"
#include "gsl/gsl_sf_coupling.h"
#include <fmt.hpp>
#include <complex>
#include <numbers>
#include <random>
#include <Eigen/Eigen>
#include "pcg/pcg_random.hpp"
#include <type_traits>


#define NO_MEMOIZE
#ifndef NO_MEMOIZE
#include <functional>
#include <map>
#include <unordered_map>
#include "tuple_hash.hpp"
#endif


namespace aef {
    [[noreturn]] inline void unreachable() {
        // Uses compiler specific extensions if possible.
        // Even if no extension is used, undefined behavior is still raised by
        // an empty function body and the noreturn attribute.
#if defined(_MSC_VER) && !defined(__clang__) // MSVC
        __assume(false);
#else // GCC, Clang
        __builtin_unreachable();
#endif
    }
    //typedef Eigen::SparseMatrix<dcomplex>;
}

extern pcg64 *pcg;
void init_rng();

template <typename T> T genrandom(T lower, T upper) {
    std::uniform_real_distribution<T> dist(lower, upper);
    return dist(*pcg);
}

inline dcomplex parity(double z) {
    using namespace std::complex_literals;
    return std::exp(1i * std::numbers::pi * z);//return std::pow(-1, dcomplex(a));
}
inline dcomplex parity(dcomplex z) {
    using namespace std::complex_literals;
    return std::exp(1i * std::numbers::pi * z);// std::pow(-1, a);
}

inline dcomplex xi(spin a, spin b) {
    return parity(a + b) * sqrt((2.0 * a + 1.0) * (2.0 * b + 1.0));
}

#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#define constexpr_sqrt sqrt
#else
#define constexpr_sqrt gcem::sqrt
#endif

constexpr double xi_prime(spin a, spin b) {
    if (std::is_constant_evaluated()) {
        return constexpr_sqrt((2 * a + 1) * (2 * b + 1));
    } else {
        return sqrt((2 * a + 1) * (2 * b + 1));
    }
}

constexpr double q_mag(spin q) {
    if (std::is_constant_evaluated()) {
        return constexpr_sqrt(q * (q + 1.0) * (2 * q + 1.0));
    } else {
        return sqrt(q * (q + 1.0) * (2 * q + 1.0));
    }
}
#include <chrono>
#include <filesystem>
namespace aef {
    std::chrono::time_point<std::chrono::system_clock> log_time_at_point(
        const char* desc,
        std::chrono::time_point<std::chrono::system_clock>& start,
        std::chrono::time_point<std::chrono::system_clock>& prev);
    bool is_aef_run_path(std::filesystem::path in);
    std::filesystem::path get_aef_run_path(std::filesystem::path in);
};
#ifndef NO_MEMOIZE
template <typename R, typename... Args>
std::function<R(Args...)> memo(R(*fn)(Args...)) {
    std::unordered_map<std::tuple<Args...>, R> table;
    return [fn, table](Args... args) mutable -> R {
        auto argt = std::make_tuple(args...);
        auto memoized = table.find(argt);
        if (memoized == table.end()) {
            auto result = fn(args...);
            table[argt] = result;
            return result;
        } else {
            return memoized->second;
        }
    };
}
#endif
/// <summary>
/// Wigner 3J coeffs
/// </summary>
/// <param name="j1"></param>
/// <param name="j2"></param>
/// <param name="j3"></param>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <param name="m3"></param>
/// <returns></returns>
static inline double w3j(double j1, double j2, double j3, double m1, double m2,
    double m3) {
    int twoj1 = (int)(2 * j1);
    int twoj2 = (int)(2 * j2);
    int twoj3 = (int)(2 * j3);

    int twom1 = (int)(2 * m1);
    int twom2 = (int)(2 * m2);
    int twom3 = (int)(2 * m3);
#ifndef NO_MEMOIZE
    static auto gsl_3j = memo(gsl_sf_coupling_3j);
    return gsl_3j(twoj1, twoj2, twoj3, twom1, twom2, twom3);
#else
    return gsl_sf_coupling_3j(twoj1, twoj2, twoj3, twom1, twom2, twom3);
#endif
}

/// <summary>
/// ClebschGordon coeffs
/// </summary>
/// <param name="j1"></param>
/// <param name="m1"></param>
/// <param name="j2"></param>
/// <param name="m2"></param>
/// <param name="J"></param>
/// <param name="M"></param>
/// <returns></returns>
static inline double cg_coeff(spin j1, spin m1, spin j2, spin m2, spin J, spin M) {
    int twoj1 = (int)(2 * j1);
    int twoj2 = (int)(2 * j2);
    int two_J = (int)(2 * J);

    int twom1 = (int)(2 * m1);
    int twom2 = (int)(2 * m2);
    int two_M = (int)(2 * M);

    double threej = gsl_sf_coupling_3j(twoj1, twoj2, two_J, twom1, twom2, -two_M);
#if 0
    double par = std::pow(-1, (twoj2 - twoj1 - two_M) / 2);
#else
    double par = std::real(parity(j2 - j1 - M));
#endif
    double degen = std::sqrt(2 * J + 1);

    return par * degen * threej;
}

static inline double w6j(double j1, double j2, double j3, double j4, double j5,
    double j6) {
    int twoj1 = (int)(2 * j1);
    int twoj2 = (int)(2 * j2);
    int twoj3 = (int)(2 * j3);

    int twoj4 = (int)(2 * j4);
    int twoj5 = (int)(2 * j5);
    int twoj6 = (int)(2 * j6);
#ifndef NO_MEMOIZE
    static auto gsl_6j = memo(gsl_sf_coupling_6j);
    return gsl_6j(twoj1, twoj2, twoj3, twoj4, twoj5, twoj6);
#else
    return gsl_sf_coupling_6j(twoj1, twoj2, twoj3, twoj4, twoj5, twoj6);
#endif
}

static inline double w9j(double j1, double j2, double j3, double j4, double j5,
    double j6, double j7, double j8, double j9) {
    int twoj1 = (int)(2 * j1);
    int twoj2 = (int)(2 * j2);
    int twoj3 = (int)(2 * j3);

    int twoj4 = (int)(2 * j4);
    int twoj5 = (int)(2 * j5);
    int twoj6 = (int)(2 * j6);

    int twoj7 = (int)(2 * j7);
    int twoj8 = (int)(2 * j8);
    int twoj9 = (int)(2 * j9);
#ifndef NO_MEMOIZE
    static auto gsl_9j = memo(gsl_sf_coupling_9j);
    return gsl_9j(twoj1, twoj2, twoj3, twoj4, twoj5, twoj6, twoj7, twoj8, twoj9);
#else
    return gsl_sf_coupling_9j(twoj1, twoj2, twoj3, twoj4, twoj5, twoj6, twoj7, twoj8, twoj9);
#endif
}

/// <summary>
/// Allows the use of std::complex<double> in fmt::format.
/// </summary>
template <> struct fmt::formatter<dcomplex> : fmt::formatter<std::string> {
    auto format(dcomplex v, format_context& ctx) const {
        return formatter<std::string>::format(
            fmt::format("({} + i*{})", std::real(v), std::imag(v)), ctx);
    }
};

/// <summary>
/// Simultaneously diagonalizes matricies A and B under the following two assumptions:
///     1) A and B are diagonalizable.
///     2) [A, B] = 0 (A and B must commute
/// 
/// Uses the algorithm provided by https://math.stackexchange.com/a/4388322: with probability 1,
/// choosing random $t\in\mathbb{R}$ and diagonalizing A+tB will provide an eigenbasis of both
/// 
/// 
/// </summary>
/// <typeparam name="Matrix"></typeparam>
/// <param name="A"></param>
/// <param name="B"></param>
/// <param name="U"></param>
/// <param name="Ut"></param>
template <class Matrix> void simultaneously_diagonalize(const Matrix& A, const Matrix& B,
    Matrix& U, Matrix& Ut) {
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
        if (agood && bgood) return;
        std::cout << fmt::format("{} FAIL FAIL FAIL", t) << std::endl;
        //DebugBreak();
    }
}

template <class Matrix, class Vector>
typename Matrix::Scalar expectation_value(const Vector &v1, const Matrix &op) {
    // for simplicity just require that scalar types are identical
    static_assert(
        std::is_same<typename Matrix::Scalar, typename Vector::Scalar>::value);
    return v1.adjoint() * op * v1;
}

template <class Matrix, class Vector>
typename Matrix::Scalar expectation_value(const Vector& v1, const Matrix& op,
    const Vector& v2) {
    // for simplicity just require that scalar types are identical
    static_assert(
        std::is_same<typename Matrix::Scalar, typename Vector::Scalar>::value);
    return v1.adjoint() * op * v2;
}

#ifdef __INTELLISENSE__
static void foo() {
    Eigen::MatrixXcd d;
    simultaneously_diagonalize(d, d, d, d);
    Eigen::VectorXcd v;
    
    expectation_value(v, d);
    expectation_value(v, d, v);
}
#endif
#endif

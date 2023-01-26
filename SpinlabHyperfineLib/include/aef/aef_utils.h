#ifndef _AEF_UTILS_H
#define _AEF_UTILS_H 1
#pragma once

#include "aef_types.h"
#include "gcem.hpp"
#include "gsl/gsl_sf_coupling.h"
#include <format>

#define NO_MEMOIZE
#ifndef NO_MEMOIZE
#include <functional>
#include <map>
#include <unordered_map>
#include "tuple_hash.hpp"
#endif

inline dcomplex parity(double a) {
    return std::pow(-1, dcomplex(a));
}
inline dcomplex parity(dcomplex a) {
    return std::pow(-1, a);
}

inline dcomplex xi(spin a, spin b) {
    return parity(a + b) * sqrt((2 * a + 1) * (2 * b + 1));
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

template <> struct std::formatter<dcomplex> : std::formatter<std::string> {
    auto format(dcomplex v, format_context& ctx) {
        return formatter<string>::format(
            std::format("({} + i*{})", std::real(v), std::imag(v)), ctx);
    }
};
#endif
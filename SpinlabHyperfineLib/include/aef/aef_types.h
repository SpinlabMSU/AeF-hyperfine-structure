#ifndef _aef_types_h
#define _aef_types_h 1

#pragma once

#include <cmath>
#include <complex>

namespace aef {
    typedef double spin;
    constexpr spin half = spin(0.5);
    typedef std::complex<double> dcomplex;
}
using aef::spin;
using aef::half;
using aef::dcomplex;

#endif
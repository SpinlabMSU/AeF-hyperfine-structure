/*
  aef/baf_constants.hpp -- constants related to the rotohyperfine structure of 138BaF
  Note: as a pure collection of constants, this file is not eligible to be copyrighted
  in the United States.  For any place/time where this is held to not be the case, 
  I place this file into the public domain following the CC0.  See 
  https://creativecommons.org/public-domain/cc0/ for details.
*/
#ifndef _BAF_CONSTANTS_HPP
#define _BAF_CONSTANTS_HPP 1
#pragma once

namespace baf_constants {
///// sourced from PRA 98, 032513 (2018)
// constants for rotational hamiltonian
constexpr double B = 6743.9586;     // MHz
constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
constexpr double gamma = 80.954;
constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

// constants for hyperfine shift
constexpr double b = 63.509; // MHz
constexpr double c = 8.224;  // MHz

// constants for stark shift
constexpr double mu_e = 3.170; // D

// constants for zeeman shift
constexpr double mu_B = 0; // magnetic moment
constexpr double mu_BN = 0; // nuclear magnetic moment of 19F
constexpr double g_S = 2; // electron spin
constexpr double g_N = 1; // nuclear
constexpr double g_r = 1; // rotational

// explicit rotational symmetry breaking term to break m_f degeneracy
// no longer needed
constexpr double e_mf_break = 0;// 1.0E-6; // MHz
};
#endif

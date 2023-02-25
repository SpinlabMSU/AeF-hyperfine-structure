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


// explicit rotational symmetry breaking term to break m_f degeneracy
// no longer needed
constexpr double e_mf_break = 0;// 1.0E-6; // MHz
};
#endif

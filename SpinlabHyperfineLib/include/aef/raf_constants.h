/*
  aef/raf_constants.hpp -- constants related to the rohyperfine structure of 225RaF
  Note: as a pure collection of constants, this file is not eligible to be copyrighted
  in the United States.  For any place/time where this is held to not be the case, 
  I place this file into the public domain following the CC0.  See 
  https://creativecommons.org/public-domain/cc0/ for details.
*/
#ifndef _RAF_CONSTANTS_HPP
#define _RAF_CONSTANTS_HPP 1
#pragma once
#include <aef/units.hpp>

namespace aef::raf_constants {
    ///// sourced from [0] PRA 98, 032513 (2018) --> these are for BaF not RaF
    //// and from [1] Nat Phys 20, 202-207 (2024)
    //// and from [2] https://arxiv.org/abs/2311.04121
    //// and from [3] Phys. Rev. A 102, 062801
    //// and from [2] https://arxiv.org/pdf/1302.5682.pdf??

    // constants for rotational hamiltonian, taken from 
    constexpr double B = 0.191985 * unit_conversion::MHz_per_inv_cm; // MHz, from [1]
    constexpr double D = 1.40E-7 * unit_conversion::MHz_per_inv_cm; // MHz, from [1]
    constexpr double gamma = 0.00585 * unit_conversion::MHz_per_inv_cm; // MHz, from [1]
    constexpr double delta = 0.111 * 1e-3; // kHz --> MHz -- not measured yet, using 138BaF value from [0]

    /// constants for Radium-225 hyperfine shift in 225RaF
    constexpr double A_para_Ra = -0.5692 * unit_conversion::MHz_per_inv_cm; // MHz, from [2]
    constexpr double A_perp_Ra = -0.5445 * unit_conversion::MHz_per_inv_cm; // MHz, from [2]
    constexpr double b_Ra = A_perp_Ra; // MHz --> used in fermi contact
    constexpr double c_Ra = A_para_Ra - A_perp_Ra;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_Ra = 0.00; // MHz --> used in nuclear spin-rotation, not measured
    // no electric quadrupole because I(225Ra) = 1/2

    /// constants for fluorine-19 hyperfine shift in 225RaF -- not currently measured,
    /// using values predicted by Phys. Rev. A 102, 062801
    constexpr double A_para_F = 109; // MHz, from [3]
    constexpr double A_perp_F = 90; // MHz, from [3]
    
    constexpr double b_F = A_perp_F; // MHz --> used in fermi contact
    constexpr double c_F = A_para_F - A_perp_F;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_F = 0.00; // MHz --> used in nuclear spin-rotation

    // no electric quadrupole because I(19F) = 1/2
    // I guess no I1 \cdot I2 term???

    // constants for stark shift
    constexpr double mu_e = 3.170; // D, from [0] for now (BaF number, not RaF)


    // explicit rotational symmetry breaking term to break m_f degeneracy
    // no longer needed
    constexpr double e_mf_break = 0;// 1.0E-6; // MHz
};

// used before 2025-08-01
namespace aef::raf_constants::old {
    ///// sourced from PRA 98, 032513 (2018) --> these are for BaF not RaF
    //// and from ??https://arxiv.org/pdf/1302.5682.pdf??
    //// and from Nat Phys 20, 202-207 (2024)
    //// and from https://arxiv.org/abs/2311.04121
    //// TODO -- look at using predictions from Phys. Rev. A 102, 062801 for 19F HFS
    // constants for rotational hamiltonian
    constexpr double B = 0.191985 * unit_conversion::MHz_per_inv_cm; // MHz
    constexpr double D = 1.40E-7 * unit_conversion::MHz_per_inv_cm; // kHz --> MHz
    constexpr double gamma = 0.00585 * unit_conversion::MHz_per_inv_cm;
    constexpr double delta = 0.111 * 1e-3; // kHz --> MHz // not measured yet, using 138BaF value

    /// constants for fluorine-19 hyperfine shift in 225RaF -- not currently measured,
    /// just using values from 138BaF for now
    constexpr double b_F = 63.509; // MHz --> used in fermi contact
    constexpr double c_F = 8.224;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_F = 0.00; // MHz --> used in nuclear spin-rotation
    // no electric quadrupole because I(19F) = 1/2

    /// constants for Radium-225 hyperfine shift in 225RaF
    constexpr double A_para_Ra = -0.5692 * unit_conversion::MHz_per_inv_cm;
    constexpr double A_perp_Ra = -0.5445 * unit_conversion::MHz_per_inv_cm;
    constexpr double b_Ra = A_perp_Ra; // MHz --> used in fermi contact
    constexpr double c_Ra = A_para_Ra - A_perp_Ra;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_Ra = 0.00; // MHz --> used in nuclear spin-rotation
    // no electric quadrupole because I(225Ra) = 1/2

    // I guess no I1 \cdot I2 term???

    // constants for stark shift
    constexpr double mu_e = 3.170; // D


    // explicit rotational symmetry breaking term to break m_f degeneracy
    // no longer needed
    constexpr double e_mf_break = 0;// 1.0E-6; // MHz
};

/// <summary>
/// This set of constants is for testing purposes only and is identical to the constants for $^{138}$BaF
/// </summary>
namespace aef::raf_constants::test {
    ///// sourced from PRA 98, 032513 (2018) --> these are for BaF not RaF
    //// and from ??https://arxiv.org/pdf/1302.5682.pdf??
    //// and from Nat Phys 20, 202-207 (2024)
    //// and from https://arxiv.org/abs/2311.04121
    // constants for rotational hamiltonian
    constexpr double B = 6743.9586;     // MHz
    constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
    constexpr double gamma = 80.954;
    constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

    /// Nucleus 1 Hyperfine structure: constants for Radium-225 hyperfine shift in pseudo-RaF
    constexpr double b_Ra = 63.509; // MHz --> used in fermi contact
    constexpr double c_Ra = 8.224;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_Ra = 0.00; // MHz --> used in nuclear spin-rotation
    // no electric quadrupole because I(225Ra) = 1/2


    /// constants for fluorine-19 hyperfine shift in 225RaF -- not currently measured,
    /// just using values from 138BaF for now
    constexpr double b_F = 0; // MHz --> used in fermi contact
    constexpr double c_F = 0;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_F = 0.00; // MHz --> used in nuclear spin-rotation
    // no electric quadrupole because I(19F) = 1/2


    // I guess no I1 \cdot I2 term???

    // constants for stark shift -- from 
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

/*
  aef/raf_constants.hpp -- constants related to the rotohyperfine structure of 225RaF
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
    ///// sourced from PRA 98, 032513 (2018) --> these are for BaF not RaF
    //// and from ??https://arxiv.org/pdf/1302.5682.pdf??
    //// and from Nat Phys 20, 202-207 (2024)
    //// and from https://arxiv.org/abs/2311.04121
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
#endif

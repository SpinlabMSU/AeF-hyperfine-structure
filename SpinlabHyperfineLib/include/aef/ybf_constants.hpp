/*
  aef/ybf_constants.hpp -- constants related to the rotohyperfine structure of ^{A}YbF
  Note: as a pure collection of constants, this file is not eligible to be copyrighted
  in the United States.  For any place/time where this is held to not be the case, 
  I place this file into the public domain following the CC0.  See 
  https://creativecommons.org/public-domain/cc0/ for details.
*/
#ifndef _YBF_CONSTANTS_HPP
#define _YBF_CONSTANTS_HPP 1
#pragma once


//////// YbFordered by nuclear spin (spin 0 first, then 1/2, then 5/2) then by abundance

namespace aef::ybf_constants {
  /// constants that are Yb-isotope independent
  // rotational constants, from PRL 74, 1554 (1995):
  constexpr double B = 6367.59; // MHz -- uncertainty is probably greater than 
  // fluorine hfs constants
  constexpr double b_F = 142.0; // MHz
  constexpr double c_F = 84.0; // MHz
  // constants for stark shift
  constexpr double mu_e = 3.91; // D
}


////////////////////////////////////////////////////////////////////
//// YbF constants for Yb isotopes with no nuclear spin         ////
//// Only eEDM measurements are sensible for these              ////
////////////////////////////////////////////////////////////////////

// 174Yb: I=0, abundance is 31.9 %
namespace aef::ybf_constants::yb_174 {
  ///// sourced from PRL 74, 1554 (1995)
  // constants for rotational hamiltonian
  constexpr double B = 6367.59;     // MHz
  constexpr double D = 5.5296 * 1e-3; // kHz --> MHz, not udpated
  constexpr double gamma = 13.31;
  constexpr double delta = -3.801 * 1e-3; // kHz --> MHz

  // constants for hyperfine shift -- 19F
  constexpr double b = 142; // MHz
  constexpr double c = 84;  // MHz
};

// 172Yb: I=0, abundance is 21.8 %
namespace aef::ybf_constants::yb_172 {
  ///// sourced from <FIND SOURCE: THIS WAS FOR 138BaF PRA 98, 032513 (2018)>
  ///// and from Phys. Chem. A 2009, 113, 28, 8038–8044 (2009)
  ///// and from PRL 74, 1554 (1995)
  // constants for rotational hamiltonian
  constexpr double B = 6367.59;     // MHz
  constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
  constexpr double gamma = 13.30;
  constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

  // constants for hyperfine shift -- 19F
  constexpr double b = 63.509; // MHz
  constexpr double c = 8.224;  // MHz

  // constants for stark shift
  constexpr double mu_e = 3.91; // D
};

// 176Yb: I=0, abundance is 12.9 %
namespace aef::ybf_constants::yb_176 {
  ///// sourced from <FIND SOURCE: THIS WAS FOR 138BaF PRA 98, 032513 (2018)>
  // constants for rotational hamiltonian
  constexpr double B = 6743.9586;     // MHz
  constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
  constexpr double gamma = 80.954;
  constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

  // constants for hyperfine shift -- 19F
  constexpr double b = 63.509; // MHz
  constexpr double c = 8.224;  // MHz

  // constants for stark shift
  constexpr double mu_e = 3.91; // D
};

// 170Yb: I=0, abundance is 3.02 %
namespace aef::ybf_constants::yb_170 {
  ///// sourced from <FIND SOURCE: THIS WAS FOR 138BaF PRA 98, 032513 (2018)>
  // constants for rotational hamiltonian
  constexpr double B = 6743.9586;     // MHz
  constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
  constexpr double gamma = 80.954;
  constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

  // constants for hyperfine shift -- 19F
  constexpr double b = 63.509; // MHz
  constexpr double c = 8.224;  // MHz

  // constants for stark shift
  constexpr double mu_e = 3.91; // D
};

////////////////////////////////////////////////////////////////////
//// YbF constants for Yb isotopes with nonzero nuclear spin    ////
//// NSM measurements are possible for for these                ////
////////////////////////////////////////////////////////////////////

// 168Yb: I=0, abundance is 0.126 %
namespace aef::ybf_constants::yb_168 {
  ///// sourced from <FIND SOURCE: THIS WAS FOR 138BaF PRA 98, 032513 (2018)>
  // constants for rotational hamiltonian
  constexpr double B = 6743.9586;     // MHz
  constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
  constexpr double gamma = 80.954;
  constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

  // constants for hyperfine shift -- 19F
  constexpr double b = 63.509; // MHz
  constexpr double c = 8.224;  // MHz

  // constants for stark shift
  constexpr double mu_e = 3.91; // D
};

// 171Yb: I=1/2, abundance is 14.2 %
namespace aef::ybf_constants::yb_171 {
    ///// sourced from PRA 98, 032513 (2018) --> these are for BaF not YbF
    //// and from ??https://arxiv.org/pdf/1302.5682.pdf?? (RaF)
    // TODO use actual YbF numbers
    // constants for rotational hamiltonian
    constexpr double B = 5689; // MHz -- approximate //6743.9586;     // MHz
    constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
    constexpr double gamma = 80.954; // MHz
    constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

    /// constants for fluorine-19 hyperfine shift in 171YbF
    constexpr double b_F = 63.509; // MHz --> used in fermi contact
    constexpr double c_F = 8.224;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_F = 0.00; // MHz --> used in nuclear spin-rotation
    // no electric quadrupole because I(19F) = 1/2

    /// constants for Yb-171 hyperfine shift in 171YbF
    constexpr double b_Yb = 63.509; // MHz --> used in fermi contact
    constexpr double c_Yb = 8.224;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_Yb = 0.00; // MHz --> used in nuclear spin-rotation
    constexpr double q_Yb = 0.00; // e*cm^2?? -- nuclear electric quadrupole

    // I guess no I1 \cdot I2 term???

    // constants for stark shift
    constexpr double mu_e = 3.91; // D
};

// 
namespace aef::ybf_constants::yb_173 {
    ///// sourced from PRA 98, 032513 (2018) --> these are for BaF not YbF
    //// and from ??https://arxiv.org/pdf/1302.5682.pdf?? (RaF)
    // TODO use actual YbF numbers
    // constants for rotational hamiltonian
    constexpr double B = 5689; // MHz -- approximate //6743.9586;     // MHz
    constexpr double D = 5.5296 * 1e-3; // kHz --> MHz
    constexpr double gamma = 13.31;
    constexpr double delta = 0.111 * 1e-3; // kHz --> MHz

    /// constants for fluorine-19 hyperfine shift in 173YbF
    constexpr double b_F = 63.509; // MHz --> used in fermi contact
    constexpr double c_F = 8.224;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_F = 0.00; // MHz --> used in nuclear spin-rotation
    // no electric quadrupole because I(19F) = 1/2

    /// constants for Yb-173 hyperfine shift in 173YbF
    constexpr double b_Yb = 63.509; // MHz --> used in fermi contact
    constexpr double c_Yb = 8.224;  // MHz --> usid in fermi contact and dipolar spin-spin
    constexpr double c_I_Yb = 0.00; // MHz --> used in nuclear spin-rotation
    constexpr double q_Yb = 0.00; // e*cm*cm?? --> nuclear electric quadrupole moment

    // I guess no I1 \cdot I2 term???

    // constants for stark shift
    constexpr double mu_e = 3.91; // D
};
#endif // _YBF_CONSTANTS_HPP
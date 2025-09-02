/*
  aef/units.hpp -- constants related to unit conversion along with a few elemental constants.
  Note: as a pure collection of constants, this file is not eligible to be copyrighted
  in the United States.  For any place/time where this is held to not be the case, 
  I place this file into the public domain following the CC0.  See 
  https://creativecommons.org/public-domain/cc0/ for details.
*/
#ifndef _UNITS_HPP
#define _UNITS_HPP 1
#pragma once

/// <summary>
/// The values of these constants are taken from NIST CODATA 2022
/// </summary>
namespace constants {
    constexpr double c = 299792458; // m/s, exact
    constexpr double h = 6.626'070'15E-34; // J*s, exact
    constexpr double k_B = 1.380'649E-23; // J/K, exact
    constexpr double mu_bohr = 13'996.244'9171; // MHz/T, +- 44 ulp (0.31 ppb) 
    constexpr double mu_nuclear = 7.622'593'2188; // MHz / T, +- 24 ulp (0.31 ppb)
    constexpr double e = 1.602176634E-19;// C, exact
};

namespace unit_conversion {
    /// <summary>
    /// Electric field conversion factor from V/m to MHz/D
    /// </summary>
    constexpr double MHz_D_per_V_m = 0.005034;
    /// <summary>
    /// Electric field conversion factor from V/cm to MHz/D
    /// </summary>
    constexpr double MHz_D_per_V_cm = 100 * MHz_D_per_V_m;
    /// <summary>
    /// The original definition of MHz_D_per_V_cm had a typo making it 10 times
    /// larger than it should have been.  This is kept here as documentation of
    /// that mistake.
    /// </summary>
    constexpr double BAD_MHz_D_per_V_cm = 1000 * MHz_D_per_V_m;
    /// <summary>
    /// Energy conversion from Rydberg to (planck's constant x) MHz 
    /// </summary>
    constexpr double MHz_per_Ry = 3.2898419603E15 / 1E6;

    /// <summary>
    /// Energy conversion from K to MHz
    /// </summary>
    constexpr double MHz_per_Kelvin = constants::k_B / constants::h * 1E-6;


    constexpr double MHz_per_inv_cm = constants::c * 1E-6 * 100;


#pragma region "Debye conversion"
    /// <summary>
    /// Conversion factor from Debye to Coulomb meters
    /// 1 Debye is defined to be 10^-18 statC*cm, which converts to 1 / c * 10^-21 C*m in SI
    /// </summary>
    constexpr double C_m_per_D = (1.0 / constants::c) * 1E-21; // C*m / D
    /// <summary>
    /// Conversion factor from Debye to e*cm
    /// </summary>
    constexpr double e_cm_per_D = C_m_per_D * 100  / constants::e; // e*cm / D
    /// <summary>
    /// Conversion factor from Debye to e*nm
    /// </summary>
    constexpr double e_nm_per_D = C_m_per_D * 1E9 / constants::e; // e*nm / D
    /// <summary>
    /// Conversion factor either from V/cm to MHz/(C*m) or from C*m to MHz / (V/cm)
    /// </summary>
    constexpr double MHz_per_V_cm_per_C_m = 1E-6 * 100/ constants::h; // MHz / (V/cm) / (C*m)
    /// <summary>
    /// Conversion factor either from V/cm to MHz/(C*m) or from C*m to MHz / (V/cm)
    /// </summary>
    constexpr double MHZ_per_V_cm_per_D = MHz_per_V_cm_per_C_m * C_m_per_D; // MHz / (V/cm) / D
}; // namespace unit_conversion

#endif //_UNITS_HPP

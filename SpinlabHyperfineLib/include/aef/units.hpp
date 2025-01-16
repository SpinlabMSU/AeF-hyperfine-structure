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

namespace constants {
    constexpr double c = 299792458; // m/s, exact
    constexpr double h = 6.626'070'15E-34; // J*s, exact
    constexpr double k_B = 1.380'649E-23; // J/K, exact
    constexpr double mu_bohr = 13'996.244'9171; // MHz/T, +- 44 ulp (0.31 ppb) 
    constexpr double mu_nuclear = 7.622'593'2188; // MHz / T, +- 24 ulp (0.31 ppb)
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
}; // namespace unit_conversion

#endif //_UNITS_HPP

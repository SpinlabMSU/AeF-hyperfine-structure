#ifndef _UNITS_HPP
#define _UNITS_HPP 1
#pragma once

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
    constexpr double MHz_per_Ry = 3.2898419603E15 / 1E6;


}; // namespace unit_conversion

#endif //_UNITS_HPP

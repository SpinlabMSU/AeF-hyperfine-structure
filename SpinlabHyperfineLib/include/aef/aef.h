/*
    aef/aef.h -- main include file for the AeF-hyperfine-structure program

    This file is part of the AeF-hyperfine-structure program. 
    
    AeF-hyperfine-structure is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version.

    AeF-hyperfine-structure is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    AeF-hyperfine-structure. If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef _AEF_H
#define _AEF_H 1
#pragma once

#ifdef _WIN32
#define _SILENCE_ALL_CXX23_DEPRECATION_WARNINGS
#define WIN32_MEAN_AND_LEAN
#include <windows.h>
#endif

#if !defined(__PRETTY_FUNCTION__) && !defined(__GNUC__)
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

namespace aef {
};

#include <numeric>
#include "aef_types.h"
#include "j_basis_vec.h"
#include "aef_utils.h"
#include "plaf_utils.hpp"
#include "baf_constants.hpp"
#include "units.hpp"
#include "HyperfineCalculator.hpp"
//#include "jf_basis_vec.h"
//#include "MolecularSystem.h"
#include "matrix_utils.h"

namespace aef {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
};
//namespace hfs_constants = baf_constants;


#endif
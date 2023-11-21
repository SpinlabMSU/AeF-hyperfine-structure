// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.
/*
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

#ifndef PCH_H
#define PCH_H

// add headers that you want to pre-compile here
#include "aef/framework.h"


#include "Eigen/Eigen"
#include "zlib.h"
#include "pcg/pcg_random.hpp"
#include "pcg/pcg_extras.hpp"

// extra external headers
#include <zstr.hpp>
#include <algorithm>
#include <iterator>
#include <gsl/gsl_sf_coupling.h>
#include <fmt.hpp>
#include <iostream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>

#endif //PCH_H

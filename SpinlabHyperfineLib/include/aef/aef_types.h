/*
    aef/aef_types.h -- contains type declarations

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
#ifndef _aef_types_h
#define _aef_types_h 1

#pragma once

#include <cmath>
#include <complex>

namespace aef {
    typedef double spin;
    constexpr spin half = spin(0.5);
    typedef std::complex<double> dcomplex;
}
using aef::spin;
using aef::half;
using aef::dcomplex;

#endif
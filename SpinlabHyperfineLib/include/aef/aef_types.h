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

#ifdef _MSC_VER
#define MSVC_ALIGN(x) __declspec(align(x))
#define ATTR_ALIGN(x)
#else
#define MSVC_ALIGN(x) 
#define ATTR_ALIGN(x) __attribute__((aligned(x)))
#endif

namespace aef {
    typedef double spin;
    constexpr spin half = spin(0.5);
    typedef MSVC_ALIGN(16) std::complex<double> dcomplex ATTR_ALIGN(16);

    /// <summary>
    /// Result codes
    /// A result code indicates success if it is non-negative
    /// </summary>
    enum class ResultCode : int32_t {
        // Successful result codes
        Success = 0, // non-detailed
        S_NOTHING_PERFORMED,
        NumSuccessfulCodes,
        /// error result codes
        OutOfMemory = -1,
        Timeout = -2,
        InvalidArgument = -3,
        IllegalState = -4,
        Unimplemented = -5,
        InternalError = -6,
        HardwareFailure = -7,
        NotAvailable = -8,
        InvalidFormat = -9,
        _LastErrorCodeP1,
        NumErrorCodes = OutOfMemory - (_LastErrorCodeP1 - 1) + 1,
        NumCodes = NumSuccessfulCodes + NumErrorCodes,
    };

    inline bool succeeded(ResultCode c) {
        return (static_cast<int32_t>(c)) >= 0;
    }

    inline bool failed(ResultCode c) {
        return !succeeded(c);
    }

}
using aef::spin;
using aef::dcomplex;

#endif

/*
    aef/operators/basis_ket.h -- contains the basis_ket concept

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
#ifndef _AEF_OPERATORS_BASIS_KET_H
#define _AEF_OPERATORS_BASIS_KET_H 1
#pragma once
#include <aef/aef_types.h>
#include <type_traits>
namespace aef::operators {
    /// <summary>
    /// The IBasisKet
    /// </summary>
    template <typename T, typename integer>
    concept IBasisKet = requires(T v, integer idx) {
        // "integer" needs to actually be an integer type
        std::integral <integer>;
        // need way of converting basis ket to integer
        {v.index()} -> std::same_as<integer>;
        // need to be able to conver 
        {T::from_index(idx)} -> std::same_as<T>;
    };
}
#endif //_AEF_OPERATORS_BASIS_KET_H
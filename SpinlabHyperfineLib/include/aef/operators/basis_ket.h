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
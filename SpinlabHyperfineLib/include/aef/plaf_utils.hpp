/*
    aef/plaf_utils.h -- contains platform-related utility functions.

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
#ifndef _AEF_PLAF_UTILS_H
#define _AEF_PLAF_UTILS_H 1
#pragma once

#include <filesystem>

/// <summary>
/// 
/// </summary>
/// <returns>The number of physical cores in the system</returns>
int get_num_cores();

bool have_read_perms(std::filesystem::path p);

inline uint32_t make_4cc(char name[4]) {
    union {
        char     in[4];
        uint32_t out;
    } conv;
    memcpy(conv.in, name, sizeof(conv.in));
    return conv.out;
}

#endif
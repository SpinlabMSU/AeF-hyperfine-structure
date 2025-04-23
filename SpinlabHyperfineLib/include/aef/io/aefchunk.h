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

#pragma once
#ifndef _AEF_IO_AEFCHUNK_H
#define _AEF_IO_AEFCHUNK_H 1

#include <stdint.h>
#include <stddef.h>

/*
 * AeFChunk is a tagged chunk-based format that I've decided to create to replace the old AeFDat format.
 * It should be substantially more flexible, make it easier for external code to generate files, and make
 * versioning/changes easier.
 */
namespace aef::chunk {
    /*
     * A traditional Four Character Code (4 ascii characters == 32 bits)
     * Note that these generally aren't zero-terminated, so don't try to use them as strings!
     */
    union fourcc {
        uint8_t cc[4];
        uint32_t ucode;

        bool operator == (fourcc other) {
            return this->ucode == other.ucode;
        }

        int operator <=> (fourcc other) {
            return ucode - other.ucode;
        }

        constexpr operator uint32_t() const {
            return std::bit_cast<uint32_t>(this->cc);
        }
    };

    /// <summary>
    /// This +
    /// </summary>
    constexpr fourcc file_magic = { .ucode = 0xAEFC'DA70 };
    /// <summary>
    /// molecular system filetype
    /// </summary>
    constexpr fourcc mol_sys = fourcc({ 'M', 'o', 'l', 55 });
    constexpr fourcc prms = fourcc({ 'p', 'r', 'm', 's' });
    constexpr fourcc oplist = fourcc({ 'O', 'p', 'L', 's' });

    /// <summary>
    /// This chunk serves as a tag marking end of file
    /// </summary>
    constexpr fourcc end0 = fourcc({ 'E', 'N', 'D', '\0' });
    /// <summary>
    /// This fourcc 
    /// </summary>
    constexpr fourcc matrix = fourcc({ 'M', 't', 'r', 'X' });

    /// <summary>
    /// Chunks with this tag store a self-adjoint operator
    /// </summary>
    constexpr fourcc self_adjoint_op = fourcc({ 0x5A, 'M', 't', 'X' });
    /// <summary>
    /// Stores the system hamiltonian and its components
    /// </summary>
    constexpr fourcc hamiltonian_chunk = fourcc({ 'H', 'a', 'M', 'C' });

    /// <summary>
    /// This struct 
    /// </summary>
    struct chunk_hdr {
        fourcc type;
        uint16_t version;
        uint16_t flags;
    };

    /// <summary>
    /// This 
    /// </summary>
    struct file_hdr {
        struct chunk_hdr hdr;
        fourcc filetype;
    };

    struct self_adjoint_operator_chunk {

    };

    struct general_matrix_chunk {
        static constexpr fourcc cc = fourcc({ 'M', 't', 'r', 'X' });
        struct chunk_hdr hdr;
        fourcc matnam;

        static constexpr fourcc nameless = { .ucode=0xFF0F'F0FF };
    };

};

namespace std {
    template <> struct hash<aef::chunk::fourcc> {
        std::size_t operator()(const aef::chunk::fourcc &fcc) const {
            return fcc.ucode;
        }
    };
};

#endif
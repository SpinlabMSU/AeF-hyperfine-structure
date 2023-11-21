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
#include "pch.h"
#include <random>
#include "aef/aef_utils.h"
#include "gsl/gsl_sf_coupling.h"
#include "pcg/pcg_extras.hpp"
#include <stdint.h>

#ifndef _WIN32
#include <unistd.h>
#endif

pcg64 *pcg;

void init_rng() {
    pcg_extras::seed_seq_from<std::random_device> sq;
    pcg = new pcg64(sq);
}



#ifdef _WIN32

#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount __popcnt
#endif
#include <tchar.h>

#elif defined(__linux__)
// we will use sysconf, so hwloc not needed
#else
#include <hwloc.h>
#endif

#define popcount __builtin_popcount


int get_num_cores() {
#ifdef _WIN32
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION pBuf = nullptr, ptr = nullptr;
    DWORD len = 0;
    bool done = false;
    
    // allocate buffer --> this is dumb
    while (!done) {
        DWORD rc = GetLogicalProcessorInformation(pBuf, &len);

        if (FALSE == rc) {
            if (GetLastError() == ERROR_INSUFFICIENT_BUFFER){
                if (pBuf)
                    free(pBuf);

                pBuf = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION)calloc(len, sizeof(char));
            } else {
                _tprintf(TEXT("\nError %d\n"), GetLastError());
                return (3);
            }
        } else {
            done = true;
        }
    }

    ptr = pBuf;

    if (!ptr) {
        std::cerr << "GetLogicalProcessorInformation returned null pointer" << std::endl;
        throw std::runtime_error("GetLogicalProcessorInformation returned null pointer in " __FILE__ ":" __PRETTY_FUNCTION__);
    }
    assert("GetLogicalProcessorInformation returned" && ptr);

    ptrdiff_t offset = 0;
    constexpr size_t sz = sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
    int numCores = 0;

    while (offset + sz <= len) {
        if (ptr->Relationship == RelationProcessorCore)
            numCores++;
        ptr++;
        offset += sz;
    }
    free(pBuf);
    return numCores;
#elif defined(__linux__)
    return sysconf(_SC_NPROCESSORS_ONLN);
#else
    int nPhysicalProcessorCount = 0;

    hwloc_topology_t sTopology;

    if (hwloc_topology_init(&sTopology) == 0 && hwloc_topology_load(sTopology) == 0) {
        nPhysicalProcessorCount = hwloc_get_nbobjs_by_type(sTopology, HWLOC_OBJ_CORE);
        hwloc_topology_destroy(sTopology);
    }
    return nPhysicalProcessorCount;
#endif
}
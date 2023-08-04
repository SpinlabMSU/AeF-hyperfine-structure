#include "pch.h"
#include <random>
#include "aef/aef_utils.h"
#include "gsl/gsl_sf_coupling.h"
#include "pcg/pcg_extras.hpp"
#include <stdint.h>

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

#else
#include <hwloc.h>
#endif

#define popcount __builtin_popcount


int get_num_cores() {
#ifdef _WIN32
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION pBuf = nullptr, ptr = nullptr;
    DWORD len;
    bool done = false;
    
    // allocate buffer --> this is dumb
    while (!done) {
        DWORD rc = GetLogicalProcessorInformation(pBuf, &len);

        if (FALSE == rc) {
            if (GetLastError() == ERROR_INSUFFICIENT_BUFFER)
            {
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
    ptrdiff_t offset = 0;
    constexpr size_t sz = sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
    int numCores = 0;

    while (offset + sz <= len) {
        if (ptr->Relationship == RelationProcessorCore)
            numCores++;
        ptr++;
        offset += sz;
    }
    delete pBuf;
    return numCores;
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
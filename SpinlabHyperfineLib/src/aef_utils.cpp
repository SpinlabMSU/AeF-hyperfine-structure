#include "pch.h"
#include <random>
#include "aef/aef_utils.h"
#include "gsl/gsl_sf_coupling.h"
#include "pcg/pcg_extras.hpp"

pcg64 *pcg;

void init_rng() {
    pcg_extras::seed_seq_from<std::random_device> sq;
    pcg = new pcg64(sq);
}
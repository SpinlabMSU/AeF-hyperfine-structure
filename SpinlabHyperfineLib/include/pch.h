// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

// add headers that you want to pre-compile here
#include "aef\framework.h"


#include "Eigen\Eigen"
#include "zlib.h"
#include "pcg/pcg_random.hpp"
#include "pcg/pcg_extras.hpp"

// extra external headers
#include <zstr.hpp>
#include <algorithm>
#include <iterator>
#include <gsl/gsl_sf_coupling.h>
#include <format>
#include <iostream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>

#endif //PCH_H

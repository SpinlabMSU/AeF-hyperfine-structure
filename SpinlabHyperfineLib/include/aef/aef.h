#ifndef _AEF_H
#define _AEF_H 1
#pragma once

#ifdef _WIN32
#define WIN32_MEAN_AND_LEAN
#include <windows.h>
#endif

#if !defined(__PRETTY_FUNCTION__) && !defined(__GNUC__)
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#include "aef_types.h"
#include "j_basis_vec.h"
#include "aef_utils.h"
#include "baf_constants.hpp"
#include "units.hpp"

namespace hfs_constants = baf_constants;


#endif
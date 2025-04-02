#ifndef _AEF_MAIN_H
#define _AEF_MAIN_H 1
#pragma once

#include <aef/aef.h>
#include <aef/debug_stream.h>
#include <aef/matrix_utils.h>
#include <aef/teestream.hpp>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fmt.hpp>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <cxxopts.hpp>
#include <aef/quantum.h>

#include <aef/MolecularSystem.h>
#include <aef/systems/BaFMolecularCalculator.h>
#include <aef/systems/RaFMolecularCalculator.h>

using namespace std::chrono;
namespace fs = std::filesystem;
namespace hfs_constants = baf_constants;
using aef::log_time_at_point;
using namespace aef::quantum;

#undef MATRIX_ELT_DEBUG

#endif // _AEF_MAIN_H
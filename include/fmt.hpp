#ifndef _FMT_HPP

#if __has_include(<format>)
#include <format>
namespace fmt = std;
#else
#include <fmt/format.h>
#include <fmt/chrono.h>
#include <fmt/ranges.h>
#include <fmt/color.h>
#endif

#endif

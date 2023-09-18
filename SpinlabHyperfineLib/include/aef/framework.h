#pragma once

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#else
static inline void DebugBreak(){}
#endif

#include <aef/aef.h>

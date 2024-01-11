/*
    aef/LogRedirector.hpp -- contains code to implement redirecting streams to the logfile.

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
#ifndef _AEF_LOG_REDIRECTOR_H
#define _AEF_LOG_REDIRECTOR_H 1
#pragma once

#include <iostream>
#include "teestream.hpp"
#include "debug_stream.h"

namespace aef {
    class LogRedirector {
    public:
        std::ostream &oLog;

        teestream::teebuf *pBuf, *pOutb, *pErrb;
        std::streambuf *pLogBuf;
        debug_stream::debug_ostream *pDebug;

        LogRedirector(std::ostream &oLog_, bool enable_debug_stream = false, bool redirect_stdio_ = true);
        ~LogRedirector();
        void touch();

    private:
        bool debug_enabled;
        bool redirect_stdio;

        std::streambuf *orig_coutb;
        std::streambuf *orig_cerrb;

        void *stdio_stdout_cookie;
        void *stdio_stdout_preserve;

        void *stdio_stderr_cookie;
        void *stdio_stderr_preserve;
    };
};
#endif

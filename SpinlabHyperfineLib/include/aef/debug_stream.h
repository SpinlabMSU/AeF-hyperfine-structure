/*
    aef/debug_stream.h -- implements a stream that outputs to the appropriate
    location for debug.  On Unix-likes this is just stderr, but on Windows one
    needs to use the OutputDebugString functions.

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
#include <streambuf>
#include <array>
#include <ostream>
#include <filesystem>

// platform-dependent
#ifdef _WIN32
#include <windows.h>
#else
#include <stdio.h>
#endif

namespace debug_stream {
    template <class charT, class traits = std::char_traits<charT>, unsigned int bufferSize = 64>
    class basic_debugbuf : public std::basic_streambuf<charT, traits> {

        using Base = std::basic_streambuf<charT, traits>;

        void writebuffer(const charT* begin, const charT* end) {
            ptrdiff_t nelem = end - begin;
            assert(nelem > 0 && nelem <= bufferSize);
            std::copy(begin, end, tbuf);
            // enforce null termination
            tbuf[bufferSize] = '\0';
#ifdef _WIN32
            if constexpr (std::is_same<traits, char>::value) {
                OutputDebugStringA(tbuf);
            } else if constexpr (std::is_same<traits, wchar_t>::value) {
                OutputDebugStringW(tbuf);
            } else {
                OutputDebugStringA((char*)tbuf);
            }
#else
            ptrdiff_t len = end - begin;
            std::cerr.write(reinterpret_cast<char*>(begin), len);
#endif
            //this->set
        }

        int overflow(int c) {
            bool rc(true);



            if (!traits::eq_int_type(traits::eof(), c)) {
                charT arr[2] = { (charT)c, (charT)0 };
#ifdef _WIN32
                if constexpr (std::is_same<traits, char>::value) {
                    OutputDebugStringA(arr);
                } else if constexpr (std::is_same<traits, wchar_t>::value) {
                    OutputDebugStringW(arr);
                } else {
                    OutputDebugStringA((char*)arr);
                }
#else
                if constexpr (std::is_same<traits, char>::value) {
                    std::cerr.write(arr, sizeof(arr));
                } else if constexpr (std::is_same<traits, wchar_t>::value) {
                    std::wcerr.put(c);
                } else {
                    std::cerr.write(arr, sizeof(arr));
                }
#endif
            }
            return rc ? traits::not_eof(c) : traits::eof();
        }

        /*
        std::streamsize xputsn(const charT* str, std::streamsize n) {
            charT *b2 =
        }
        */
        int sync() {
            // always synced
            return 0;
        }

    private:
        std::array<charT, 64> _buffer;
        std::array <charT, bufferSize + 1> tbuf;
    public:
        basic_debugbuf() : _buffer{}, tbuf{} {
            //Base::setp(_buffer.begin(), _buffer.end());
        }
    };

    typedef basic_debugbuf<char> debugbuf;

    template <class charT, class traits = std::char_traits<charT>> class basic_debug_ostream
        : private virtual basic_debugbuf<charT, traits>
        , public std::basic_ostream<charT, traits> {
    public:
        basic_debug_ostream() : std::basic_ostream<charT, traits>(this) {
            this->init(this);
        }
    };

    typedef basic_debug_ostream<char> debug_ostream;
};

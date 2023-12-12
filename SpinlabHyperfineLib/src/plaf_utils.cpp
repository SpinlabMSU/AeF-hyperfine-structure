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
#include <iostream>
using namespace std::chrono;
namespace fs = std::filesystem;

#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#else
#include <io.h>
#define access _access
#define F_OK 0
#define W_OK (1<<1)
#define R_OK (1<<2)
#endif

bool have_read_perms(std::filesystem::path p) {
    const char *path = (char*)(p.c_str());
    errno = 0;
    if (!access(path, F_OK)) {
        assert(errno == EACCES || errno == ENOENT || errno == ENOTDIR);
        // file/dir doesn't exist
        return false;
    }

    if (!access(path, R_OK)) {
        assert(errno == EACCES || errno == ENOENT || errno == ENOTDIR);
        // no read access
        return false;
    }

    if (fs::is_directory(p)) {
        // can't open directories with normal open on most modern OSes.
        return true;
    }
    // simplest test is to just try to open the file
    // putting this in a seperate function allows for 
    return (bool)std::ifstream(p);
}
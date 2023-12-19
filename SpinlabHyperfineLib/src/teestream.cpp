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
#include "aef/teestream.hpp"

namespace teestream {
    namespace fs = std::filesystem;

    teebuf *tee_cout(fs::path& opath)
    {
        std::ofstream* oLog = new std::ofstream(opath, std::ios::out | std::ios::trunc);

        std::streambuf* obuf = std::cout.rdbuf();
        teebuf* pBuff = new teebuf(obuf, oLog->rdbuf());
        std::cout.rdbuf(pBuff);//std::cout.set_rdbuf(pBuff);
        return pBuff;
    }

    teebuf *tee_cerr(fs::path& opath)
    {
        std::ofstream* oLog = new std::ofstream(opath, std::ios::out | std::ios::trunc);

        std::streambuf* obuf = std::cerr.rdbuf();
        teebuf* pBuff = new teebuf(obuf, oLog->rdbuf());
        std::cerr.rdbuf(pBuff);
        return pBuff;
    }

    teebuf *tee_cout_cerr(fs::path& opath) {
        std::ofstream *oLog = new std::ofstream(opath, std::ios::out | std::ios::trunc);
        
        std::streambuf *obuf = std::cout.rdbuf();
        teebuf *pBuff = new teebuf(obuf, oLog->rdbuf());
        std::cout.rdbuf(pBuff);

        tee_cerr(opath);

        return pBuff;
    }
};

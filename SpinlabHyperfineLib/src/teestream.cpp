#include "pch.h"
#include "aef/teestream.hpp"

namespace tee {
    namespace fs = std::filesystem;

    teebuf *tee_cout(fs::path& opath)
    {
        std::ofstream* oLog = new std::ofstream(opath, std::ios::out | std::ios::trunc);

        std::streambuf* obuf = std::cout.rdbuf();
        teebuf* pBuff = new teebuf(obuf, oLog->rdbuf());
        std::cout.set_rdbuf(pBuff);
        return pBuff;
    }

    teebuf *tee_cerr(fs::path& opath)
    {
        std::ofstream* oLog = new std::ofstream(opath, std::ios::out | std::ios::trunc);

        std::streambuf* obuf = std::cerr.rdbuf();
        teebuf* pBuff = new teebuf(obuf, oLog->rdbuf());
        std::cerr.set_rdbuf(pBuff);
    }

    teebuf *tee_cout_cerr(fs::path& opath) {
        std::ofstream *oLog = new std::ofstream(opath, std::ios::out | std::ios::trunc);
        
        std::streambuf *obuf = std::cout.rdbuf();
        teebuf *pBuff = new teebuf(obuf, oLog->rdbuf());
        std::cout.set_rdbuf(pBuff);

        tee_cerr(opath);

        return pBuff;
    }
};
#include "pch.h"
#include "aef/LogRedirector.h"

aef::LogRedirector::LogRedirector(std::ostream &oLog_, bool enable_debug_stream, bool redirect_stdio_) : oLog(oLog_) {
    debug_enabled = enable_debug_stream;
    redirect_stdio = redirect_stdio_;

    orig_coutb = std::cout.rdbuf();
    orig_cerrb = std::cerr.rdbuf();

    pDebug = new debug_stream::debug_ostream;

    // test if logging to debug window enabled
    if (enable_debug_stream) {
        // if enabled, tee to both the log file and the debug window
        pBuf = new tee::teebuf(oLog.rdbuf(), pDebug->rdbuf());
        pLogBuf = pBuf;
    } else {
        // otherwise just output to the log file
        pBuf = nullptr;
        pLogBuf = oLog.rdbuf();
    }
    pOutb = new tee::teebuf(pLogBuf, orig_coutb);
    std::cout.rdbuf(pOutb);

    pErrb = new tee::teebuf(pLogBuf, orig_cerrb);
    std::cerr.rdbuf(pErrb);
}

aef::LogRedirector::~LogRedirector() {
    oLog.flush();

    if (redirect_stdio) {
        fflush(stdout);
        fflush(stderr);
        // todo redirection
    }

    // redirect
    std::cout.rdbuf(orig_coutb);
    std::cout.rdbuf(orig_cerrb);

    // delete teebufs
    delete pOutb; pOutb = nullptr;
    delete pErrb; pErrb = nullptr;
    // don't delete pLogBuf becuase we don't own oLog

}

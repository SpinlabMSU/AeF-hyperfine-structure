#include "pch.h"
#include "aef/LogRedirector.h"

#ifdef __GLIBC__
#include <stdio.h>
#include <iostream>
namespace {
    struct aef_stdio_cookie {
        cookie_io_functions_t funcs;
        std::ostream *os;
        FILE **redirect_stream;
        FILE *stream;
        FILE *saved_stream;
        bool closed;

        aef_stdio_cookie(std::ostream *os_, FILE **redirect_stream_): 
            os(os_), redirect_stream(redirect_stream_) {
            funcs.read = reinterpret_cast<cookie_read_function_t*>(&aef_stdio_cookie::read);
            funcs.write = reinterpret_cast<cookie_write_function_t*>(&aef_stdio_cookie::write);
            funcs.seek = nullptr;
            funcs.close = reinterpret_cast<cookie_close_function_t*>(&aef_stdio_cookie::close);

            stream = fopencookie(this, "w", funcs);
            saved_stream = *redirect_stream;
            *redirect_stream = stream;

            closed = false;
        }

        ~aef_stdio_cookie(){
            if (!closed) {
                aef_stdio_cookie::close(this);
            }
            //memset((void*)this, 0, sizeof(*this));
        }
        static ssize_t read(aef_stdio_cookie *c,char *buf, size_t size) {
            return -1;
        }

        static ssize_t write(aef_stdio_cookie *c, const char *buf, size_t size) {
            c->os->write(buf, size);
            return size;
        }

        static int close(aef_stdio_cookie *c){
            *(c->redirect_stream) = c->saved_stream;
            c->closed = true;
            return 0;
        }
    };
};
#endif

aef::LogRedirector::LogRedirector(std::ostream &oLog_, bool enable_debug_stream, bool redirect_stdio_) : oLog(oLog_) {
    debug_enabled = enable_debug_stream;
    redirect_stdio = redirect_stdio_;

    orig_coutb = std::cout.rdbuf();
    orig_cerrb = std::cerr.rdbuf();

    pDebug = new debug_stream::debug_ostream;

    // test if logging to debug window enabled
    if (enable_debug_stream) {
        // if enabled, tee to both the log file and the debug window
        pBuf = new teestream::teebuf(oLog.rdbuf(), pDebug->rdbuf());
        pLogBuf = pBuf;
    } else {
        // otherwise just output to the log file
        pBuf = nullptr;
        pLogBuf = oLog.rdbuf();
    }
    pOutb = new teestream::teebuf(pLogBuf, orig_coutb);
    std::cout.rdbuf(pOutb);

    pErrb = new teestream::teebuf(pLogBuf, orig_cerrb);
    std::cerr.rdbuf(pErrb);
#ifdef __GLIBC__
    aef_stdio_cookie *stdout_cook = new aef_stdio_cookie(&std::cout, &stdout);
    stdio_stdout_cookie = (void*) stdout_cook;
    aef_stdio_cookie *stderr_cook = new aef_stdio_cookie(&std::cerr, &stderr);
    stdio_stderr_cookie = (void*) stderr_cook;
#endif
}

aef::LogRedirector::~LogRedirector() {
    oLog.flush();
// temporarily prevent crashes
#ifdef AEF_ACTUALLY_CLOSE_THIS
    if (redirect_stdio) {
        fflush(stdout);
        fflush(stderr);
        // todo redirection
#ifdef __GLIBC__
    aef_stdio_cookie *stdout_cook = (aef_stdio_cookie *)stdio_stdout_cookie;
    delete stdout_cook; stdio_stdout_cookie = nullptr;
    aef_stdio_cookie *stderr_cook = (aef_stdio_cookie *)stdio_stderr_cookie;
    delete stderr_cook; stdio_stderr_cookie = nullptr;
#endif
    }

    // redirect
    std::cout.rdbuf(orig_coutb);
    std::cerr.rdbuf(orig_cerrb);

    // delete teebufs
    delete pOutb; pOutb = nullptr;
    delete pErrb; pErrb = nullptr;
    // don't delete pLogBuf becuase we don't own oLog
#endif
}


#if defined(__GNUC__) && (defined(__x86_64__) || defined(__I386__)) 
#include <cpuid.h>
#endif
void aef::LogRedirector::touch() {
#if defined(__GNUC__)
    #if defined(__I386__) || defined(__x86_64__)
    // serializing inst
    //asm volatile("SERIALIZE": : :"memory"); // crashes with SIGILL
    asm volatile("": : :"memory");
    // use cpuid instead
    union {
      char text[16];
      uint32_t reg[4];
    } mfr;
    explicit_bzero(&mfr, sizeof(mfr));
    uint32_t max_leaf = 0;
    // leaf 0: max_leaf in eax, ManufacturerID in {ebx,edx,ecx} 
    __cpuid(0, max_leaf, mfr.reg[0], mfr.reg[2], mfr.reg[1]);
    std::cout << fmt::format("CPU manufacturer {}, max leaf id {}", mfr.text, max_leaf) << std::endl;
    #elif defined(__aarch64__)
    // 
    asm volatile("dmb ish\ndsb ish\nisb sy": : :"memory");
    #endif
    __sync_synchronize();
#elif defined(_MSC_VER)
    int cpuInfo[4];
    __cpuid(cpuInfo, 0);
    char cpu_mfr[sizeof(cpuInfo) + 1];
    SecureZeroMemory(cpu_mfr, sizeof(cpu_mfr));
    memcpy(cpu_mfr, cpuInfo, sizeof(cpuInfo) - sizeof(int));
    std::cout << "Printing manufacturer info to prevent spurious optimization :" << cpu_mfr << std::endl;
    MemoryBarrier();
#endif
}
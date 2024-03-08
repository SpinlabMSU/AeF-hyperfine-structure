#include "pch.h"
#include "aef/matrix_utils.h"
#include "aef/backends/EigenMatrixBackend.h"
#ifndef DONT_USE_CUDA
#include "aef/backends/CudaMatrixBackend.h"
#include "cuda_utils.h"
#endif
#include <atomic>

namespace aef::matrix {
    namespace {
        IMatrixOpBackend* backend = nullptr;
        std::atomic_int initcount = 0;
    };
    ResultCode init(BackendType hint, int argc, char** argv) {
        if (initcount++ > 0) {
            assert(backend != nullptr);
            return ResultCode::S_NOTHING_PERFORMED;
        }
        ResultCode res = ResultCode::NotAvailable;
        // select backend
#ifndef DONT_USE_CUDA
        // try cuda first
        CudaMatrixBackend* cu = new CudaMatrixBackend();
        res = cu->init(argc, argv);
        // success --> set and return
        if (static_cast<int32_t>(res) > 0) {
            backend = cu;
            return res;
        }
        // 
#endif
        backend = get_fallback_backend();
        return backend->init(argc, argv);
    }
    IMatrixOpBackend* aef::matrix::get_backend() {
        return backend;
    }
    IMatrixOpBackend* get_fallback_backend() {
        static EigenMatrixBackend back;
        return &back;
    }
    ResultCode shutdown() {
        // not initialized
        if (backend == nullptr) {
            return ResultCode::S_NOTHING_PERFORMED;
        }
        // init count > 0
        if (--initcount > 0) {
            return ResultCode::S_NOTHING_PERFORMED;
        }
        ResultCode res = backend->shutdown();
        if (backend != get_fallback_backend()) {
            delete backend;
        }
        backend = nullptr;
        return res;
    }
}
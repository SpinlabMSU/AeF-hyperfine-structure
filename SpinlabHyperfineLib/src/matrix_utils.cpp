#include "pch.h"
#include "aef/matrix_utils.h"
#include "aef/backends/EigenMatrixBackend.h"
#ifndef DONT_USE_CUDA
#include "aef/backends/CudaMatrixBackend.h"
#include "cuda_utils.h"
#endif
#include "aef/backends/OneMklMatrixBackend.h"
#include <atomic>

namespace aef::matrix {
    namespace {
        IMatrixOpBackend* backend = nullptr;
        std::atomic_int initcount = 0;
    };

    bool is_init() {
        return backend != nullptr;
    }

    ResultCode init(BackendType hint, int argc, char** argv) {
        if (initcount++ > 0) {
            assert(backend != nullptr);
            return ResultCode::S_NOTHING_PERFORMED;
        }
        ResultCode res = ResultCode::NotAvailable;

        switch (hint) {
        default:
        case BackendType::NvidiaCuda:
#ifndef DONT_USE_CUDA
            ;
            backend = new CudaMatrixBackend();
            res = backend->init(argc, argv);

            if (succeeded(res)) {
                return res;
            }
            delete backend;
            backend = nullptr;
            [[fallthrough]];
#else
            // cuda disabled --> can't use
            [[fallthrough]];
#endif
        case BackendType::IntelOneAPI:;
            backend = new OneMklMatrixBackend();
            res = backend->init(argc, argv);

            if (succeeded(res)) {
                return res;
            }
            delete backend;
            backend = nullptr;
        case BackendType::EigenCPU:
            backend = get_fallback_backend();
            return backend->init(argc, argv);
        }
        // 
        aef::unreachable();
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

    ResultCode commutes(bool& out, Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd* workspace, double prec) {
        bool need_alloc = !workspace;
        if (need_alloc) {
            workspace = new Eigen::MatrixXcd();
            workspace->resizeLike(A);
        }

        ResultCode res = aef::matrix::commutator(A, B, *workspace);

        if (failed(res)) {
            goto end;
        }
        out = workspace->isZero(prec);
        end:
        if (need_alloc) {
            delete workspace;
        }
        return res;
    }

    bool commutes(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd *workspace, double prec) {
        bool result = false;
        ResultCode res = commutes(result, A, B, workspace, prec);
        return result;
    }
}
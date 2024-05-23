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
        std::cout << "[aef::matrix] Selecting backend" << std::endl;
        switch (hint) {
        default:
        case BackendType::NvidiaCuda:
#ifndef DONT_USE_CUDA
            std::cout << "[aef::matrix] Trying Cuda backend" << std::endl;
            backend = new CudaMatrixBackend();
            res = backend->init(argc, argv);

            if (succeeded(res)) {
                return res;
            }
            std::cout << "[aef::matrix] Cuda backend failed, trying next" << std::endl;
            delete backend;
            backend = nullptr;
            [[fallthrough]];
#else
            // cuda disabled --> can't use
            [[fallthrough]];
#endif
        case BackendType::IntelOneAPI:
            std::cout << "[aef::matrix] Trying Intel OneAPI backend" << std::endl;;
            backend = new OneMklMatrixBackend();
            res = backend->init(argc, argv);

            if (succeeded(res)) {
                return res;
            }
            std::cout << "[aef::matrix] Intel OneAPI backend failed, trying next" << std::endl;
            delete backend;
            backend = nullptr;
            [[fallthrough]];
        case BackendType::EigenCPU:
            std::cout << "[aef::matrix] Using EigenCPU backend.  This cannot fail." << std::endl;
            backend = get_fallback_backend();
            return backend->init(argc, argv);
        }
        // 
        std::cout << "[aef::matrix] Unreachable hit, this should never happen." << std::endl;
        aef::unreachable();
    }
    IMatrixOpBackend* get_backend() {
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
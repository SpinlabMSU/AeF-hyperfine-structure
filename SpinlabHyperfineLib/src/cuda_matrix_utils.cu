//#include "pch.h"

#define _SILENCE_ALL_CXX23_DEPRECATION_WARNINGS
#define _AEF_WILL_USE_CUDA_HEADERS
#include "aef/matrix_utils.h"

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <cusolverDn.h>
#include "cuda_utils.h"
#include "cusolver_utils.h"

namespace aef {
    // cuda and cusolver props
    namespace {
        bool init = false;
        // selected device ID
        int devID = -1;
        // stream
        cudaStream_t cu_stream;
        cusolverDnHandle_t cu_handle = 0;
        cudaDeviceProp deviceProps;
        // device matrix pointer --> used both for input hermitian matrix and for evec output
        cuDoubleComplex* d_A;
        // device eigenvalues pointer --> note: this is real bec
        double *d_W;
        // device info ptr
        int *d_info = nullptr;

        //
        int saved_n = -1;
    };


    bool init_cuda(int argc, const char **argv) {
        if (init) {
            return true;
        }
        devID = cuda::findCudaDevice(argc, argv);
        checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
        checkCudaErrors(cudaSetDevice(devID));
        checkCudaErrors(cusolverDnCreate(&cu_handle));
        checkCudaErrors(cudaStreamCreateWithFlags(&cu_stream, cudaStreamNonBlocking));
        checkCudaErrors(cusolverDnSetStream(cu_handle, cu_stream));
        saved_n = -1;
        init = true;
        return true;
    }

    bool shutdown_cuda() {
        if (!init) {
            return true;
        }
        cudaStreamSynchronize(cu_stream);
        if (d_A) {
            checkCudaErrors(cudaFreeAsync(d_A, cu_stream));
            d_A = nullptr;
        }

        if (d_W) {
            checkCudaErrors(cudaFreeAsync(d_W, cu_stream));
            d_W = nullptr;
        }

        if (d_info) {
            checkCudaErrors(cudaFreeAsync(d_info, cu_stream));
            d_info = nullptr;
        }

        checkCudaErrors(cudaStreamSynchronize(cu_stream));
        checkCudaErrors(cusolverDnDestroy(cu_handle));
        cu_handle = nullptr;
        checkCudaErrors(cudaStreamDestroy(cu_stream));
        cu_stream = nullptr;
        checkCudaErrors(cudaDeviceReset());
        init = false;
        return true;
    }

    bool is_cuda_initialized() {
        return init;
    }

    void mat_init(cudaStream_t stream) {
        auto status = cusolverDnCreate(&cu_handle);
        cu_stream = stream;
        CUSOLVER_CHECK(cusolverDnSetStream(cu_handle, stream));
        init = true;
    }

    void cuda_resize(int n) {
        assert(n >= 0);
        if (n == saved_n) {
            return;
        }
        const size_t szA = sizeof(cuDoubleComplex) * n * n;
        const size_t szW = sizeof(double) * n;

        if (d_A) {
            checkCudaErrors(cudaFreeAsync(d_A, cu_stream));
            d_A = nullptr;
        }

        if (d_W) {
            checkCudaErrors(cudaFreeAsync(d_W, cu_stream));
            d_W = nullptr;
        }

        if (d_info) {
            checkCudaErrors(cudaFreeAsync(d_info, cu_stream));
            d_info = nullptr;
        }

        if (n == 0) {
            // don't bother allocating zero-sized arrays
            saved_n = n;
            return;
        }
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_A), szA, cu_stream));
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_W), szW, cu_stream));
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_info), sizeof(int), cu_stream));

        checkCudaErrors(cudaStreamSynchronize(cu_stream));
        saved_n = n;
    }

    void log_dev_props_info(std::ostream& out) {
        //deviceProps.
    }


    void diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) {
        const int rows = (int)mat.rows();

        if (rows <= 0) {
            // don't do work on a zero-sized matrix
            return;
        }

        std::cout << "[Cuda-based diagonalizer] Diagonalize called" << std::endl;

        if (rows > saved_n || saved_n <= 0) {
            cuda_resize(rows);
        }

        const size_t mat_size = sizeof(cuDoubleComplex) * mat.size();
        const size_t ws_size = sizeof(cuDoubleComplex) * rows;
        const cuDoubleComplex *data = reinterpret_cast<cuDoubleComplex*>(mat.data());
        cuDoubleComplex *pW = reinterpret_cast<cuDoubleComplex*>(evals.data());
        cuDoubleComplex *pV = reinterpret_cast<cuDoubleComplex*>(evecs.data());
        int info = 0;

        // upload to GPU
        checkCudaErrors(cudaMemcpyAsync(d_A, data, mat_size, cudaMemcpyHostToDevice, cu_stream));
        std::cout << "[Cuda-based diagonalizer] data uploaded to gpu" << std::endl;
        // allocate workspace: first query how large it needs to be, then allocate
        const auto jobz = CUSOLVER_EIG_MODE_VECTOR;
        const auto uplo = CUBLAS_FILL_MODE_UPPER;
        int lwork;
        checkCudaErrors(cusolverDnZheevd_bufferSize(cu_handle, jobz, uplo, rows, d_A, rows, d_W, &lwork));
        std::cout << "[Cuda-based diagonalizer] zheev work size will be " << lwork * sizeof(cuDoubleComplex)<< " bytes" << std::endl;
        cuDoubleComplex* d_Work;
        checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_Work), lwork * sizeof(cuDoubleComplex)));
        std::cout << "[Cuda-based diagonalizer] allocated work space on gpu" << std::endl;
        
        // call cusolvers ZHEEV, then copy data back to CPU ram
        auto status = (cusolverDnZheevd(cu_handle, jobz, uplo, rows, d_A, rows, d_W, d_Work, lwork, d_info));
        std::cout << "[Cuda-based diagonalizer] queued zheev execution" << std::endl;
        checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost, cu_stream));
        std::cout << "[Cuda-based diagonalizer] scheduled zheev info output to be copied back to host" << std::endl;
        if (info != 0 || status != CUSOLVER_STATUS_SUCCESS) {
            // errcode
            cudaStreamSynchronize(cu_stream);
            std::cout << "cuSOLVER ZHEEV execution failed in " __FILE__ " at line # " << 164 << " info is " << info << std::endl;
            checkCudaErrors(status);
        }

        checkCudaErrors(cudaMemcpyAsync(pV, d_A, mat_size, cudaMemcpyDeviceToHost, cu_stream));
        checkCudaErrors(cudaMemcpyAsync(pW, d_W, ws_size , cudaMemcpyDeviceToHost, cu_stream));
        std::cout << "[Cuda-based diagonalizer] scheduled for data to be copied back to host" << std::endl;
        // wait for all operations to complete
        checkCudaErrors(cudaStreamSynchronize(cu_stream));
        std::cout << "[Cuda-based diagonalizer] diagonalizaion has completed execution" << std::endl;
        
        //
        if (info != 0 || status != CUSOLVER_STATUS_SUCCESS) {
            // errcode
            std::cout << "cuSOLVER ZHEEV execution failed in " __FILE__ " at line # " << 164 << " info is " << info << std::endl;
            checkCudaErrors(status);
        }
    }
}
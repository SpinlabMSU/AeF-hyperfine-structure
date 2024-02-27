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
        cublasHandle_t hCublas = nullptr;
        cudaDeviceProp deviceProps;
        // device matrix pointer --> used both for input hermitian matrix and for evec output
        cuDoubleComplex* d_A; // main matrix ptr
        cuDoubleComplex* d_U; // used for unitary
        // device eigenvalues pointer --> note: this is real bec
        double *d_W;
        cuDoubleComplex* d_V;
        int lwork = 0;
        cuDoubleComplex* d_Work = nullptr; // also used as temproary
        // host eigenvalues --> need this because 
        std::vector<double> h_W;
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
        checkCudaErrors(cublasCreate(&hCublas));
        checkCudaErrors(cublasSetStream(hCublas, cu_stream));
        saved_n = -1;
        init = true;
        return true;
    }

    static inline void ensure_work_capacity(size_t num_elts) {
        if (num_elts <= lwork) {
            return;
        }
        std::cout << "[aef::cuda_utils] need to reallocate work buffer, "
            "old size was " << lwork * sizeof(cuDoubleComplex) << " bytes, new size will be "
            << num_elts * sizeof(cuDoubleComplex) << " bytes" << std::endl;
        checkCudaErrors(cudaFree(d_Work));
        checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_Work), num_elts* sizeof(cuDoubleComplex)));
        lwork = num_elts;
        std::cout << "[aef::cuda_utils] allocated new work space on gpu" << std::endl;
    }

    bool shutdown_cuda() {
        if (!init) {
            return true;
        }
        cudaStreamSynchronize(cu_stream);

        #define freePtr(d_X) if (d_X) { checkCudaErrors(cudaFreeAsync(d_X, cu_stream)); d_X = nullptr;}
        freePtr(d_A);
        freePtr(d_U);
        freePtr(d_W);
        freePtr(d_V);
        freePtr(d_info);
        freePtr(d_Work);
        #undef freePtr

        checkCudaErrors(cudaStreamSynchronize(cu_stream));
        checkCudaErrors(cusolverDnDestroy(cu_handle));
        cu_handle = nullptr;
        checkCudaErrors(cublasDestroy(hCublas));
        hCublas = nullptr;
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
        cu_stream = stream;
        checkCudaErrors(cusolverDnCreate(&cu_handle));
        CUSOLVER_CHECK(cusolverDnSetStream(cu_handle, stream));
        checkCudaErrors(cublasCreate(&hCublas));
        checkCudaErrors(cublasSetStream(hCublas, cu_stream));
        saved_n = -1;
        init = true;
    }

    void cuda_resize(int n) {
        assert(n >= 0);
        if (n == saved_n) {
            return;
        }

        assert(init);

        std::cout << "[Cuda matrix backend] Resizing from " << saved_n << " rows to " << n << " rows." << std::endl;
        const size_t szV = sizeof(cuDoubleComplex) * n;
        const size_t szA = sizeof(cuDoubleComplex) * n * n;
        const size_t szW = sizeof(double) * n;
        size_t szTotal = szV + 2 * szA + szW;
        std::cout << "[Cuda matrix backend] Estimated initial allocation size is " << szTotal << "bytes = " << szTotal / (1 << 20) << "MiB" << std::endl;

        if (d_A) {
            checkCudaErrors(cudaFreeAsync(d_A, cu_stream));
            d_A = nullptr;
        }

        if (d_W) {
            checkCudaErrors(cudaFreeAsync(d_W, cu_stream));
            d_W = nullptr;
        }

        if (d_V) {
            checkCudaErrors(cudaFreeAsync(d_V, cu_stream));
            d_V = nullptr;
        }

        if (d_info) {
            checkCudaErrors(cudaFreeAsync(d_info, cu_stream));
            d_info = nullptr;
        }

        if (d_U) {
            checkCudaErrors(cudaFreeAsync(d_U, cu_stream));
        }
        if (n == 0) {
            // don't bother allocating zero-sized arrays
            h_W.resize(0);
            saved_n = n;
            return;
        }
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_V), szV, cu_stream));
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_A), szA, cu_stream));
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_U), szA, cu_stream));
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_W), szW, cu_stream));
        CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_info), sizeof(int), cu_stream));

        // allocate host-side eigenvalue buffer
        h_W.reserve(n);
        // pre-allocate workspace: first query how large it needs to be, then allocate
        const auto jobz = CUSOLVER_EIG_MODE_VECTOR;
        const auto uplo = CUBLAS_FILL_MODE_UPPER;
        checkCudaErrors(cusolverDnZheevd_bufferSize(cu_handle, jobz, uplo, n, d_A, n, d_W, &lwork));
        const size_t szWork = lwork * sizeof(cuDoubleComplex);
        std::cout << "[Cuda matrix backend] zheev work size will be " << szWork << " bytes" << std::endl;
        szTotal += szWork;
        std::cout << "[Cuda matrix backend] Estimated total allocation size is " << szTotal << "bytes = " << szTotal / (1 << 20) << "MiB" << std::endl;
        checkCudaErrors(cudaMallocAsync(reinterpret_cast<void**>(&d_Work), lwork * sizeof(cuDoubleComplex), cu_stream));

        checkCudaErrors(cudaStreamSynchronize(cu_stream));
        saved_n = n;
        std::cout << "[Cuda matrix backend] Resizing to size " << n << " complete" << std::endl;
    }

    void log_dev_props_info(std::ostream& out) {
        //deviceProps.
    }

    dcomplex cuda_expectation_value(Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::MatrixXcd& v2) {
        const int rows = A.rows();
        const size_t As_size = sizeof(cuDoubleComplex) * A.size();
        const size_t vs_size = sizeof(double) * A.rows();

        constexpr auto uplo = CUBLAS_FILL_MODE_UPPER;
        constexpr cuDoubleComplex one = { 1.0, 0.0 };
        constexpr cuDoubleComplex zero = { 0.0, 0.0 };

        checkCudaErrors(cudaMemcpyAsync(d_V, v1.data(), vs_size, cudaMemcpyHostToDevice, cu_stream));
        checkCudaErrors(cudaMemcpyAsync(d_A,  A.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
        //checkCudaErrors(cublasZsymv(hCublas, uplo, &one, d_A, rows, d_V, 1, &zero, nullptr))

        return dcomplex();
    }

    dcomplex cuda_expectation_value(Eigen::VectorXcd& v1, Eigen::MatrixXcd& A) {
        const size_t mat_size = sizeof(cuDoubleComplex) * A.size();
        const size_t ws_size = sizeof(double) * A.rows();


        return dcomplex();
    }

    void cuda_expectation_values(Eigen::MatrixXcd& U, Eigen::MatrixXcd& A, Eigen::MatrixXcd &out) {
        const size_t As_size = sizeof(cuDoubleComplex) * A.size();
        const int rows = A.rows();

        constexpr auto nop = CUBLAS_OP_N;
        constexpr auto oph = CUBLAS_OP_C;
        constexpr cuDoubleComplex one = { 1.0, 0.0 };
        constexpr cuDoubleComplex zero = { 0.0, 0.0 };
        // ensure work buffer is large enough
        ensure_work_capacity(A.size());
        // memcpy
        checkCudaErrors(cudaMemcpyAsync(d_A, A.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
        checkCudaErrors(cudaMemcpyAsync(d_U, U.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
        checkCudaErrors(cublasZgemm(hCublas, nop, nop, rows, rows, rows, &one, d_A, rows, d_U, rows, &zero, d_Work, rows));
        checkCudaErrors(cublasZgemm(hCublas, oph, nop, rows, rows, rows, &one, d_U, rows, d_Work, rows, &zero, d_A, rows));
        checkCudaErrors(cudaMemcpyAsync(out.data(), d_A, As_size, cudaMemcpyDeviceToHost, cu_stream));
        checkCudaErrors(cudaStreamSynchronize(cu_stream));
    }


    void diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) {
        const int rows = (int)mat.rows();

        if (rows <= 0) {
            // don't do work on a zero-sized matrix
            return;
        }

        std::cout << "[Cuda-based diagonalizer] Diagonalize called" << std::endl;

        if (rows > saved_n || saved_n <= 0) {
            std::cout << "[Cuda-based diagonalizer] Automatically resizing from " << saved_n << " rows to " << rows << " rows." << std::endl;
            cuda_resize(rows);
        }

        const size_t mat_size = sizeof(cuDoubleComplex) * mat.size();
        const size_t ws_size = sizeof(double) * rows;
        const cuDoubleComplex *data = reinterpret_cast<cuDoubleComplex*>(mat.data());
        double* pW = h_W.data();//reinterpret_cast<cuDoubleComplex*>(evals.data());
        cuDoubleComplex *pV = reinterpret_cast<cuDoubleComplex*>(evecs.data());
        int info = 0;

        // upload to GPU
        checkCudaErrors(cudaMemcpyAsync(d_A, data, mat_size, cudaMemcpyHostToDevice, cu_stream));
        std::cout << "[Cuda-based diagonalizer] data uploaded to gpu" << std::endl;
        // check workspace buffer size is large enough
        const auto jobz = CUSOLVER_EIG_MODE_VECTOR;
        const auto uplo = CUBLAS_FILL_MODE_UPPER;
        int job_lwork = 0;
        checkCudaErrors(cusolverDnZheevd_bufferSize(cu_handle, jobz, uplo, rows, d_A, rows, d_W, &job_lwork));
        // reallocate if necessary
        if (job_lwork > lwork) {
            std::cout << "[Cuda-based diagonalizer] need to reallocate zheev work space, "
                "old size was " << lwork * sizeof(cuDoubleComplex) << " bytes, new size will be " 
                << job_lwork * sizeof(cuDoubleComplex) << " bytes" << std::endl;
            checkCudaErrors(cudaFree(d_Work));
            checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_Work), job_lwork * sizeof(cuDoubleComplex)));
            lwork = job_lwork;
            std::cout << "[Cuda-based diagonalizer] allocated new work space on gpu" << std::endl;
        } else {
            std::cout << "[Cuda-based diagonalizer] using pre-allocated workspace of " << lwork * sizeof(cuDoubleComplex) << " bytes." << std::endl;
        }
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

        // copy from host memory to eigenvalue vector
        // this is necessary because evals is a dcomplex vector
        // but CUDA outputs a real double vector.
        //std::copy(h_W.begin(), h_W.end(), evals.data());
        for (int i = 0; i < evals.size(); i++) {
            evals(i) = h_W[i];
        }
        std::cout << "[Cuda-based diagonalizer] fixed up eigenvalue vector" << std::endl;

        //
        if (info != 0 || status != CUSOLVER_STATUS_SUCCESS) {
            // errcode
            std::cout << "cuSOLVER ZHEEV execution failed in " __FILE__ " at line # " << 164 << " info is " << info << std::endl;
            checkCudaErrors(status);
        }
    }
}
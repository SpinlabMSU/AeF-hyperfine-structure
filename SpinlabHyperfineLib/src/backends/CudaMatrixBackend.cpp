#define _SILENCE_ALL_CXX23_DEPRECATION_WARNINGS
#define _AEF_WILL_USE_CUDA_HEADERS
// not defining this causes a bunch of errors about cusparse using deprecated types to declare deprecated APIs
// when using MSVC.  
#define DISABLE_CUSPARSE_DEPRECATED
#include "pch.h"
#include "aef/backends/CudaMatrixBackend.h"

#include "../cuda_utils.h"
#include "../cusolver_utils.h"

using aef::matrix::ResultCode;


aef::matrix::CudaMatrixBackend::CudaMatrixBackend():
    deviceProps() {
    cu_handle = nullptr;
    cu_stream = nullptr;
    hCublas = nullptr;
    d_A = nullptr;
    d_U = nullptr;
    d_W = nullptr;
    d_V = nullptr;
    _init = false;
    devID = -1;
}

aef::matrix::CudaMatrixBackend::~CudaMatrixBackend() {
}

ResultCode aef::matrix::CudaMatrixBackend::init(int argc, char** argv) {
    if (_init) {
        return ResultCode::S_NOTHING_PERFORMED;
    }
    devID = cuda::findCudaDevice(argc, const_cast<const char**>(argv));

    if (devID < 0) {
        return ResultCode::NotAvailable;
    }

    checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
    checkCudaErrors(cudaSetDevice(devID));
    checkCudaErrors(cusolverDnCreate(&cu_handle));
    checkCudaErrors(cudaStreamCreateWithFlags(&cu_stream, cudaStreamNonBlocking));
    checkCudaErrors(cusolverDnSetStream(cu_handle, cu_stream));
    checkCudaErrors(cublasCreate(&hCublas));
    checkCudaErrors(cublasSetStream(hCublas, cu_stream));
    saved_n = -1;

    int val;
    checkCudaErrors(cudaDeviceGetAttribute(&val, cudaDevAttrMemoryPoolsSupported, devID));
    std::cout << "Cuda attribute memory pool supported is " << val << " for device " << devID << std::endl;

    _init = true;
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::ensureWorkCapacity(size_t num_elts) {
    if (!_init) {
        return ResultCode::IllegalState;
    }
    if (num_elts <= lWork) {
        return ResultCode::S_NOTHING_PERFORMED;
    }
    std::cout << "[aef::cuda_utils] need to reallocate work buffer, "
        "old size was " << lWork * sizeof(cuDoubleComplex) << " bytes, new size will be "
        << num_elts * sizeof(cuDoubleComplex) << " bytes" << std::endl;
    checkCudaErrors(cudaFree(d_Work));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_Work), num_elts * sizeof(cuDoubleComplex)));
    lWork = num_elts;
    std::cout << "[aef::cuda_utils] allocated new work space on gpu" << std::endl;
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::shutdown() {
    if (!_init) {
        return ResultCode::S_NOTHING_PERFORMED;
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
    _init = false;
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::set_max_size(int nMaxDim) {
    if (!_init) {
        return ResultCode::IllegalState;
    }
    int n = nMaxDim;
    assert(n >= 0);
    if (n == saved_n) {
        return ResultCode::S_NOTHING_PERFORMED;
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

#ifdef NO_MATRIX_ALLOCATION_HACK
    if (d_U) {
        checkCudaErrors(cudaFreeAsync(d_U, cu_stream));
    }
#endif
    if (n == 0) {
        // don't bother allocating zero-sized arrays
        h_W.resize(0);
        saved_n = n;
        return ResultCode::Success;
    }
    CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_V), szV, cu_stream));
    CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_A), szA, cu_stream));
#ifdef NO_MATRIX_ALLOCATION_HACK
    CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_U), szA, cu_stream));
#endif
    CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_W), szW, cu_stream));
    CUDA_CHECK(cudaMallocAsync(reinterpret_cast<void**>(&d_info), sizeof(int), cu_stream));

    // allocate host-side eigenvalue buffer
    h_W.reserve(n);
    // pre-allocate workspace: first query how large it needs to be, then allocate
    const auto jobz = CUSOLVER_EIG_MODE_VECTOR;
    const auto uplo = CUBLAS_FILL_MODE_UPPER;
    checkCudaErrors(cusolverDnZheevd_bufferSize(cu_handle, jobz, uplo, n, d_A, n, d_W, &lWork));
#ifdef NO_MATRIX_ALLOCATION_HACK
    const size_t szWork = lWork * sizeof(cuDoubleComplex);
    std::cout << "[Cuda matrix backend] zheev work size will be " << szWork << " bytes, " << lWork << " elements." << std::endl;
    szTotal += szWork;
#else
    // Matrix allocation hack: 
    const size_t szWork = std::max(lWork * sizeof(cuDoubleComplex), 2 * szA);
    const bool bOverride = (szWork > lWork * sizeof(cuDoubleComplex));
    const char* strOverride = bOverride ? "overrided by matrix hack" : "not overrided by matrix hack";
    std::cout << "[Cuda matrix backend] zheev work size will be " << szWork << " bytes, " << strOverride << ", " <<
        szA << " bytes re-used for d_U." << std::endl;
    if (bOverride) {
        const int lWorkNew = (szWork + sizeof(cuDoubleComplex) - 1) / sizeof(cuDoubleComplex);
        std::cout << "[Cuda matrix backend]" << "zheev minimum-required lWork was " << lWork << " elements, will now be " << lWorkNew << " elements" << std::endl;
        lWork = lWorkNew;
    } else {
        std::cout << "[Cuda matrix backend] lWork is " << lWork << " elements" << std::endl;
    }
    szTotal += (szWork - szA);
#endif
    std::cout << "[Cuda matrix backend] Estimated total allocation size is " << szTotal << "bytes = " << szTotal / (1 << 20) << "MiB" << std::endl;
    checkCudaErrors(cudaMallocAsync(reinterpret_cast<void**>(&d_Work), szWork, cu_stream));
#ifndef NO_MATRIX_ALLOCATION_HACK
    d_U = d_Work + szA;
#endif
    checkCudaErrors(cudaStreamSynchronize(cu_stream));
    saved_n = n;
    std::cout << "[Cuda matrix backend] Resizing to size " << n << " complete" << std::endl;
    
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::multiply(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
    if (!_init) {
        return ResultCode::IllegalState;
    }
    const size_t As_size = sizeof(cuDoubleComplex) * A.size();
    const int rows = A.rows();

    constexpr auto nop = CUBLAS_OP_N;
    constexpr auto oph = CUBLAS_OP_C;
    constexpr cuDoubleComplex one = { 1.0, 0.0 };
    constexpr cuDoubleComplex zero = { 0.0, 0.0 };
    // ensure work buffer is large enough
    ensureWorkCapacity(A.size());
    // memcpy
    checkCudaErrors(cudaMemcpyAsync(d_A, A.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cudaMemcpyAsync(d_U, B.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cublasZgemm(hCublas, nop, nop, rows, rows, rows, &one, d_A, rows, d_U, rows, &zero, d_Work, rows));
    checkCudaErrors(cudaMemcpyAsync(out.data(), d_Work, As_size, cudaMemcpyDeviceToHost, cu_stream));
    checkCudaErrors(cudaStreamSynchronize(cu_stream));
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::commutator(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
    if (!_init) {
        return ResultCode::IllegalState;
    }
    const size_t As_size = sizeof(cuDoubleComplex) * A.size();
    const int rows = A.rows();

    constexpr auto nop = CUBLAS_OP_N;
    constexpr auto oph = CUBLAS_OP_C;
    constexpr cuDoubleComplex one = { 1.0, 0.0 };
    constexpr cuDoubleComplex zero = { 0.0, 0.0 };
    constexpr cuDoubleComplex mone = { 0.0, -1.0 };
    // ensure work buffer is large enough
    ensureWorkCapacity(A.size());
    // memcpy
    checkCudaErrors(cudaMemcpyAsync(d_A, A.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cudaMemcpyAsync(d_U, B.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cublasZgemm(hCublas, nop, nop, rows, rows, rows, &one, d_A, rows, d_U, rows, &zero, d_Work, rows));
    checkCudaErrors(cublasZgemm(hCublas, oph, nop, rows, rows, rows, &one, d_U, rows, d_A, rows, &mone, d_Work, rows));
    checkCudaErrors(cudaMemcpyAsync(out.data(), d_Work, As_size, cudaMemcpyDeviceToHost, cu_stream));
    checkCudaErrors(cudaStreamSynchronize(cu_stream));
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::group_action(Eigen::MatrixXcd& out, Eigen::MatrixXcd& U, Eigen::MatrixXcd& A) {
    if (!_init) {
        return ResultCode::IllegalState;
    }
    const size_t As_size = sizeof(cuDoubleComplex) * A.size();
    const int rows = A.rows();

    constexpr auto nop = CUBLAS_OP_N;
    constexpr auto oph = CUBLAS_OP_C;
    constexpr cuDoubleComplex one = { 1.0, 0.0 };
    constexpr cuDoubleComplex zero = { 0.0, 0.0 };
    // ensure work buffer is large enough
    ensureWorkCapacity(A.size());
    // memcpy
    checkCudaErrors(cudaMemcpyAsync(d_A, A.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cudaMemcpyAsync(d_U, U.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cublasZgemm(hCublas, nop, nop, rows, rows, rows, &one, d_A, rows, d_U, rows, &zero, d_Work, rows));
    checkCudaErrors(cublasZgemm(hCublas, oph, nop, rows, rows, rows, &one, d_U, rows, d_Work, rows, &zero, d_A, rows));
    checkCudaErrors(cudaMemcpyAsync(out.data(), d_A, As_size, cudaMemcpyDeviceToHost, cu_stream));
    checkCudaErrors(cudaStreamSynchronize(cu_stream));
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A) {
    if (!_init) {
        return ResultCode::IllegalState;
    }

    if (v1.rows() != A.cols() || A.rows() != A.cols()) {
        // can't multiply 
        return ResultCode::InvalidArgument;
    }

    const size_t As_size = sizeof(cuDoubleComplex) * A.size();
    const size_t vs_size = sizeof(cuDoubleComplex) * v1.size();
    const int dim = A.rows();

    constexpr auto nop = CUBLAS_OP_N;
    constexpr auto oph = CUBLAS_OP_C;
    constexpr cuDoubleComplex one = { 1.0, 0.0 };
    constexpr cuDoubleComplex zero = { 0.0, 0.0 };
    cuDoubleComplex result;

    checkCudaErrors(cudaMemcpyAsync(d_A, A.data() , As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cudaMemcpyAsync(d_V, v1.data(), vs_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cublasZgemv(hCublas, nop, dim, dim, &one, d_A, dim, d_V, 1, &zero, d_Work, 1));
    checkCudaErrors(cublasZdotc(hCublas, dim, d_V, 1, d_Work, 1, &result));
    checkCudaErrors(cudaStreamSynchronize(cu_stream));
    out = { result.x, result.y };
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::matrix_element(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::VectorXcd& v2) {
    if (!_init) {
        return ResultCode::IllegalState;
    }

    if (v1.rows() != A.cols() || A.rows() != v2.rows()) {
        // can't multiply 
        return ResultCode::InvalidArgument;
    }

    const size_t As_size = sizeof(cuDoubleComplex) * A.size();
    const size_t vs_size = sizeof(cuDoubleComplex) * v1.size();
    const int dim = A.rows();

    constexpr auto nop = CUBLAS_OP_N;
    constexpr auto oph = CUBLAS_OP_C;
    constexpr cuDoubleComplex one = { 1.0, 0.0 };
    constexpr cuDoubleComplex zero = { 0.0, 0.0 };
    cuDoubleComplex result;

    checkCudaErrors(cudaMemcpyAsync(d_A, A.data(), As_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cudaMemcpyAsync(d_V, v1.data(), vs_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cudaMemcpyAsync(d_U, v2.data(), vs_size, cudaMemcpyHostToDevice, cu_stream));
    checkCudaErrors(cublasZgemv(hCublas, nop, dim, dim, &one, d_A, dim, d_V, 1, &zero, d_Work, 1));
    checkCudaErrors(cublasZdotc(hCublas, dim, d_U, 1, d_Work, 1, &result));
    checkCudaErrors(cudaStreamSynchronize(cu_stream));
    out = { result.x, result.y };
    return ResultCode::Success;
}

ResultCode aef::matrix::CudaMatrixBackend::diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) {
    if (!_init) {
        return ResultCode::IllegalState;
    }
    const int rows = (int)mat.rows();

    if (rows <= 0) {
        // don't do work on a zero-sized matrix
        return ResultCode::S_NOTHING_PERFORMED;
    }

    std::cout << "[aef::matrix::CudaMatrixBackend] Diagonalize called" << std::endl;

    if (rows > saved_n || saved_n <= 0) {
        std::cout << "[aef::matrix::CudaMatrixBackend] Automatically resizing from " << saved_n << " rows to " << rows << " rows." << std::endl;
        cuda_resize(rows);
    }

    const size_t mat_size = sizeof(cuDoubleComplex) * mat.size();
    const size_t ws_size = sizeof(double) * rows;
    const cuDoubleComplex* data = reinterpret_cast<cuDoubleComplex*>(mat.data());
    double* pW = h_W.data();
    cuDoubleComplex* pV = reinterpret_cast<cuDoubleComplex*>(evecs.data());
    int info = 0;

    // upload to GPU
    checkCudaErrors(cudaMemcpyAsync(d_A, data, mat_size, cudaMemcpyHostToDevice, cu_stream));
    std::cout << "[aef::matrix::CudaMatrixBackend] data uploaded to gpu" << std::endl;
    // check workspace buffer size is large enough
    const auto jobz = CUSOLVER_EIG_MODE_VECTOR;
    const auto uplo = CUBLAS_FILL_MODE_UPPER;
    int job_lwork = 0;
    checkCudaErrors(cusolverDnZheevd_bufferSize(cu_handle, jobz, uplo, rows, d_A, rows, d_W, &job_lwork));
    // reallocate if necessary
    if (job_lwork > lWork) {
        std::cout << "[aef::matrix::CudaMatrixBackend] need to reallocate zheev work space, "
            "old size was " << lWork * sizeof(cuDoubleComplex) << " bytes, new size will be "
            << job_lwork * sizeof(cuDoubleComplex) << " bytes" << std::endl;
        checkCudaErrors(cudaFree(d_Work));
        checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_Work), job_lwork * sizeof(cuDoubleComplex)));
        lWork = job_lwork;
        std::cout << "[aef::matrix::CudaMatrixBackend] allocated new work space on gpu" << std::endl;
    } else {
        std::cout << "[aef::matrix::CudaMatrixBackend] using pre-allocated workspace of " << lWork * sizeof(cuDoubleComplex) << " bytes." << std::endl;
    }
    // call cusolvers ZHEEV, then copy data back to CPU ram
    auto status = (cusolverDnZheevd(cu_handle, jobz, uplo, rows, d_A, rows, d_W, d_Work, lWork, d_info));
    std::cout << "[aef::matrix::CudaMatrixBackend] queued zheev execution" << std::endl;
    checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost, cu_stream));
    std::cout << "[aef::matrix::CudaMatrixBackend] scheduled zheev info output to be copied back to host" << std::endl;
    if (info != 0 || status != CUSOLVER_STATUS_SUCCESS) {
        // errcode
        cudaStreamSynchronize(cu_stream);
        std::cout << "cuSOLVER ZHEEV execution failed in " __FILE__ " at line # " << 164 << " info is " << info << std::endl;
        checkCudaErrors(status);
    }

    checkCudaErrors(cudaMemcpyAsync(pV, d_A, mat_size, cudaMemcpyDeviceToHost, cu_stream));
    checkCudaErrors(cudaMemcpyAsync(pW, d_W, ws_size, cudaMemcpyDeviceToHost, cu_stream));
    std::cout << "[aef::matrix::CudaMatrixBackend] scheduled for data to be copied back to host" << std::endl;
    // wait for all operations to complete
    checkCudaErrors(cudaStreamSynchronize(cu_stream));
    std::cout << "[aef::matrix::CudaMatrixBackend] diagonalizaion has completed execution" << std::endl;

    // copy from host memory to eigenvalue vector
    // this is necessary because evals is a dcomplex vector
    // but CUDA outputs a real double vector.
    //std::copy(h_W.begin(), h_W.end(), evals.data());
    for (int i = 0; i < evals.size(); i++) {
        evals(i) = h_W[i];
    }
    std::cout << "[aef::matrix::CudaMatrixBackend] fixed up eigenvalue vector" << std::endl;

    //
    if (info != 0 || status != CUSOLVER_STATUS_SUCCESS) {
        // errcode
        std::cout << "cuSOLVER ZHEEV execution failed in " __FILE__ " at line # " << 164 << " info is " << info << std::endl;
        checkCudaErrors(status);
    }
    return ResultCode::Success;
}
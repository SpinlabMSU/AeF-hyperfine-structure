#pragma once
#include "aef/MatrixOpBackend.h"

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <cusolverDn.h>

namespace aef::matrix {
    class CudaMatrixBackend : public IMatrixOpBackend {
        // has initialize been called
        bool _init = false;
        // selected device ID
        int devID = -1;
        // Cuda handles
        cudaStream_t cu_stream;
        cusolverDnHandle_t cu_handle;
        cublasHandle_t hCublas;
        cudaDeviceProp deviceProps;
        // device matrix pointer --> used both for input hermitian matrix and for evec output
        cuDoubleComplex* d_A; // main matrix ptr
        cuDoubleComplex* d_U; // used for unitary
        // device eigenvalues pointer --> note: this is real because CUDA outputs that way
        double* d_W;
        cuDoubleComplex* d_V;
        int lWork = 0;
        cuDoubleComplex* d_Work = nullptr; // also used as temproary
        // host eigenvalues --> need this because 
        std::vector<double> h_W;
        // device info ptr
        int* d_info = nullptr;

        //
        int saved_n = -1;

        ResultCode ensureWorkCapacity(size_t nElements);

    public:
        // Note: constructor does not 
        CudaMatrixBackend();
        ~CudaMatrixBackend();
        /// <summary>
        /// Initialize the backend with optional arguments, specified the same
        /// </summary>
        /// <param name="argc">Argument count</param>
        /// <param name="argv">Argument vector, can be null if argc == 0</param>
        /// <returns>Result code</returns>
        virtual ResultCode init(int argc, char** argv);
        /// <summary>
        /// Close out the backend
        /// </summary>
        /// <returns>Result code</returns>
        virtual ResultCode shutdown();

        /// <summary>
        /// Allocate enough space for 
        /// </summary>
        /// <param name="nMaxDim">Maximum matrix dimension</param>
        /// <returns>Result code</returns>
        virtual ResultCode set_max_size(int nMaxDim);

        /// <summary>
        /// Multiplies two matricies
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="out"></param>
        /// <returns>Result code</returns>
        virtual ResultCode multiply(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out);

        /// <summary>
        /// Computes the commutator [A, B]
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="out"></param>
        /// <returns></returns>
        virtual ResultCode commutator(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out);

        /// <summary>
        /// Computes the "group action" UAU^{-1} of 
        /// </summary>
        /// <param name="U">A unitary matrix</param>
        /// <param name="A">A general matrix</param>
        /// <param name="out">The output matrix</param>
        /// <returns>Result code</returns>
        virtual ResultCode group_action(Eigen::MatrixXcd& out, Eigen::MatrixXcd& U, Eigen::MatrixXcd& A);
        /// <summary>
        /// Calculates the expectation value of operator A in state v1
        /// </summary>
        /// <param name="out">expectation value</param>
        /// <param name="v1"></param>
        /// <param name="A"></param>
        /// <returns></returns>
        virtual ResultCode expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A);
        /// <summary>
        /// Calculates the matrix element <v
        /// </summary>
        /// <param name="out"></param>
        /// <param name="v1"></param>
        /// <param name="A"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        virtual ResultCode matrix_element(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::VectorXcd& v2);

        /// <summary>
        /// Diagonalizes complex Hermitian matricies (ZHEEV)
        /// </summary>
        /// <param name="mat">The matrix to diagonalize, must be hermitian</param>
        /// <param name="evals"></param>
        /// <param name="evecs"></param>
        virtual ResultCode diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs);
    };

}
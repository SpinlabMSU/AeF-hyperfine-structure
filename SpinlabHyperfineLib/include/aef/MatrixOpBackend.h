#pragma once
#ifndef _AEF_MATRIX_OP_BACKEND_H
#define _AEF_MATRIX_OP_BACKEND_H 1
#include "aef_types.h"
#ifndef __NVCC__
#include "aef_utils.h"
#endif
#include "gcem.hpp"
#include <complex>
#include <numbers>
#include <Eigen/Eigen>

namespace aef::matrix {
    using aef::ResultCode;
    /// <summary>
    /// Interface describing what operations must be implemented by a backend
    /// </summary>
    class IMatrixOpBackend {
    public:
        /// <summary>
        /// Initialize the backend with optional arguments, specified the same
        /// </summary>
        /// <param name="argc">Argument count</param>
        /// <param name="argv">Argument vector, can be null if argc == 0</param>
        /// <returns>Result code</returns>
        virtual ResultCode init(int argc, char** argv) = 0;
        /// <summary>
        /// Close out the backend
        /// </summary>
        /// <returns>Result code</returns>
        virtual ResultCode shutdown() = 0;

        /// <summary>
        /// Allocate enough space for 
        /// </summary>
        /// <param name="nMaxDim">Maximum matrix dimension</param>
        /// <returns>Result code</returns>
        virtual ResultCode set_max_size(int nMaxDim) = 0;

        /// <summary>
        /// Multiplies two matricies
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="out"></param>
        /// <returns>Result code</returns>
        virtual ResultCode multiply(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) = 0;

        /// <summary>
        /// Computes the commutator [A, B]
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="out"></param>
        /// <returns></returns>
        virtual ResultCode commutator(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) = 0;

        /// <summary>
        /// Computes the "group action" UAU^{-1} of 
        /// </summary>
        /// <param name="U">A unitary matrix</param>
        /// <param name="A">A general matrix</param>
        /// <param name="out">The output matrix</param>
        /// <returns>Result code</returns>
        virtual ResultCode group_action(Eigen::MatrixXcd& out, Eigen::MatrixXcd& U, Eigen::MatrixXcd& A) = 0;
        /// <summary>
        /// Calculates the expectation value of operator A in state v1
        /// </summary>
        /// <param name="out">expectation value</param>
        /// <param name="v1"></param>
        /// <param name="A"></param>
        /// <returns></returns>
        virtual ResultCode expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A) = 0;
        /// <summary>
        /// Calculates the matrix element <v
        /// </summary>
        /// <param name="out"></param>
        /// <param name="v1"></param>
        /// <param name="A"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        virtual ResultCode matrix_element(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::VectorXcd& v2) = 0;

        /// <summary>
        /// Diagonalizes complex Hermitian matricies (ZHEEV)
        /// </summary>
        /// <param name="mat">The matrix to diagonalize, must be hermitian</param>
        /// <param name="evals"></param>
        /// <param name="evecs"></param>
        virtual ResultCode diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) = 0;
    };

    /// <summary>
    /// Potential Backends -- right now only
    /// </summary>
    enum class BackendType:int32_t {
        Invalid = 0,
        // default CPU-based backend
        EigenCPU = 1,
        // CUDA -- will be the only single-device-vendor backend supported
        NvidiaCuda = 2,
        // Intel OneAPI: probable next backend implementation since it works on all main vendors
        IntelOneAPI = 3,
        OpenCL = 4,
        VKCompute = 5,
        // Will not implement
        // AMDRocmSockem -- no matter how good the pun sounds, this only works on _some_ AMD GPUs --> no
        // AppleMetal -- I don't have any devices supporting this, and it only works on Apple's devices --> no
        // 
        NumBackends
    };

    // TODO add backend chooser stuff
};

#endif //_AEF_MATRIX_OP_BACKEND_H
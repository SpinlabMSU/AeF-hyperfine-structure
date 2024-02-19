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


    enum class StatusCode : int32_t {
        Success = 0,
        OutOfMemory = 1,
        Timeout = 2,
        InvalidArgument = 3,
        IllegalState = 4,
        Unimplemented = 5,
        InternalError = 6,
        HardwareFailure = 7,
        NotAvailable = 8,
        NumCodes
    };

    /// <summary>
    /// Interface describing what operations must be implemented by a backend
    /// </summary>
    class IMatrixOpBackend {
        virtual StatusCode init(int argc, char** argv) = 0;
        virtual StatusCode shutdown() = 0;

        virtual StatusCode resize(int nmax) = 0;

        virtual StatusCode multiply(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) = 0;
        /// <summary>
        /// Computes UAU^{-1}
        /// </summary>
        /// <param name="U"></param>
        /// <param name="A"></param>
        /// <param name="out"></param>
        /// <returns></returns>
        virtual StatusCode expectation_values(Eigen::MatrixXcd& out, Eigen::MatrixXcd& U, Eigen::MatrixXcd& A) = 0;
        virtual StatusCode expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A) = 0;
        virtual StatusCode expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::VectorXcd& v2) = 0;

        /// <summary>
        /// Diagonalizes complex Hermitian matricies (ZHEEV)
        /// </summary>
        /// <param name="mat">The matrix to diagonalize, must be hermitian</param>
        /// <param name="evals"></param>
        /// <param name="evecs"></param>
        virtual StatusCode diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) = 0;
    };

    /// <summary>
    /// Potential Backends -- right now only
    /// </summary>
    enum class BackendTypes:int32_t {
        Invalid = 0,
        // default CPU-based backend
        EigenCPU = 1,
        // CUDA -- will be the only single-device-vendor backend supported
        NvidiaCuda = 2,
        // Intel OneAPI: probable next backend implementation since it works on all main vendors
        IntelOneAPI = 3,
        OpenCL = 4,
        // Will not implement
        // AMDRocmSockem -- no matter how good the pun sounds, this only works on _some_ AMD GPUs --> no
        // AppleMetal -- I don't have any devices supporting this, and it only works on Apple's devices --> no
        // 
        NumBackends
    };

    // TODO add backend chooser stuff
};

#endif //_AEF_MATRIX_OP_BACKEND_H
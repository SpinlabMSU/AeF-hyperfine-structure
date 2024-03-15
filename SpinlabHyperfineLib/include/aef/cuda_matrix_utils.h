#pragma once

#include "aef_types.h"
#ifndef __NVCC__
#include "aef_utils.h"
#endif
#include "gcem.hpp"
#include <complex>
#include <numbers>
#include <Eigen/Eigen>

namespace aef {
    /// <summary>
    /// Initializes cuda
    /// </summary>
    /// <param name="argc">Argument count</param>
    /// <param name="argv">Argument list pointer</param>
    /// <returns>True if initializing cuda succeded, false otherwise (note: as currently implemented, actually just crashes)</returns>
    bool init_cuda(int argc, const char** argv);
    /// <summary>
    /// Shuts down cuda
    /// </summary>
    /// <returns></returns>
    bool shutdown_cuda();

    bool is_cuda_initialized();

    /// <summary>
    /// Sets the size of the GPU-side CUDA buffers
    /// </summary>
    /// <param name="n"></param>
    void cuda_resize(int n);

    void log_dev_props_info(std::ostream& out);

    dcomplex cuda_expectation_value(Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::MatrixXcd& v2);
    dcomplex cuda_expectation_value(Eigen::VectorXcd& v1, Eigen::MatrixXcd& A);
    void cuda_expectation_values(Eigen::MatrixXcd& U, Eigen::MatrixXcd& A, Eigen::MatrixXcd& out);

    /// <summary>
    /// 
    /// </summary>
    /// <param name="mat"></param>
    /// <param name="evals"></param>
    /// <param name="evecs"></param>
    void diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs);
}
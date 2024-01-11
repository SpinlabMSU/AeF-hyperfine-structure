/* Copyright (c) 2022, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 ////////////////////////////////////////////////////////////////////////////////
 // These are CUDA Helper functions for initialization and error checking

#ifndef COMMON_HELPER_CUDA_H_
#define COMMON_HELPER_CUDA_H_

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <cusparse_v2.h>
#include <cusolverDn.h>
#include <cusolver_common.h>
#include <curand.h>

#include "helper_string.h"
#ifndef MAX
#define MAX(a, b) (a > b ? a : b)
#endif

// This will output the proper CUDA error strings in the event
// that a CUDA host call returns an error
#define checkCudaErrors(val) cuda::check((val), #val, __FILE__, __LINE__)

// This will output the proper error string when calling cudaGetLastError
#define getLastCudaError(msg) cuda::__getLastCudaError(msg, __FILE__, __LINE__)

// This will only print the proper error string when calling cudaGetLastError
// but not exit program incase error detected.
#define printLastCudaError(msg) cuda::__printLastCudaError(msg, __FILE__, __LINE__)

/// <summary>
/// This namespace provides declarations related to CUDA
/// </summary>
namespace cuda {
    const char* cudaGetErrorEnum(cudaError_t error);
    const char* cudaGetErrorEnum(CUresult error);
    const char* cudaGetErrorEnum(cublasStatus_t error);
    const char* cudaGetErrorEnum(cufftResult error);
    const char* cudaGetErrorEnum(cusparseStatus_t error);
    const char* cudaGetErrorEnum(cusolverStatus_t error);
    const char* cudaGetErrorEnum(curandStatus_t error);
    // const char* cudaGetErrorEnum(nvjpegStatus_t error);
    // const char* cudaGetErrorEnum(NppStatus error);

    // Float To Int conversion
    inline int ftoi(float value) {
        return (value >= 0 ? static_cast<int>(value + 0.5)
            : static_cast<int>(value - 0.5));
    }

    // General GPU Device CUDA Initialization
    int gpuDeviceInit(int devID);
    // This function returns the best GPU (with maximum GFLOPS)
    int gpuGetMaxGflopsDeviceId();
    // Initialization code to find the best CUDA Device
    int findCudaDevice(int argc, const char** argv);
    /// <summary>
    /// find integrated gpu supporting cuda, if any exist
    /// </summary>
    int findIntegratedGPU();
    // General check for CUDA GPU SM Capabilities
    bool checkCudaCapabilities(int major_version, int minor_version);

    /// <summary>
    /// Error checking
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="result"></param>
    /// <param name="func"></param>
    /// <param name="file"></param>
    /// <param name="line"></param>
    template <typename T>
    void check(T result, char const* const func, const char* const file,
        int const line) {
        if (result) {
            fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
                static_cast<unsigned int>(result), cudaGetErrorEnum(result), func);
            exit(EXIT_FAILURE);
        }
    }

    // mostly-internal helpers
    void __getLastCudaError(const char* errorMessage, const char* file,
        const int line);
    void __printLastCudaError(const char* errorMessage, const char* file,
        const int line);
    const char* _ConvertSMVer2ArchName(int major, int minor);
    int _ConvertSMVer2Cores(int major, int minor);
};
#endif  // COMMON_HELPER_CUDA_H_
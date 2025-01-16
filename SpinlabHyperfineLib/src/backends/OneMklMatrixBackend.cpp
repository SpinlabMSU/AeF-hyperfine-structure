#include "pch.h"
#include "aef/backends/OneMklMatrixBackend.h"
#ifdef HAVE_MKL
#include "oneapi/mkl.hpp"
#else
// just enough to compile for right now 
namespace sycl{
    void *gpu_selector(){
        return nullptr;
    }
    struct device{
        device(void *v){
            (void)v;
        }
    };
    struct queue{
        queue(device& d){
            (void)d;
        }
    };
}
#endif

using aef::matrix::ResultCode;

struct aef::matrix::OneMKLdata {
    bool init;
#ifdef _USE_ONEAPI
    sycl::device* device;
    sycl::queue *queue;
#else
    int* device;
    int* queue;
#endif
    dcomplex* dev_A;
    dcomplex* dev_U;

    OneMKLdata() {
        device = nullptr;
        queue = nullptr;
        sizeof(*device);
        sizeof(*queue);
    }
};

aef::matrix::OneMklMatrixBackend::OneMklMatrixBackend() {
    ptr = new OneMKLdata;
}

aef::matrix::OneMklMatrixBackend::~OneMklMatrixBackend() {
    delete ptr;
    ptr = nullptr;
}

ResultCode aef::matrix::OneMklMatrixBackend::init(int argc, char** argv) {
    if (ptr->init) {
        return ResultCode::S_NOTHING_PERFORMED;
    }
#ifdef _USE_ONEAPI
    ptr->device = new sycl::device(sycl::gpu_selector());
    ptr->queue = new sycl::queue(*(ptr->device));
#else
    // need fewer ifdefs by doing this
    ptr->device = new int;
    ptr->queue = new int;
#endif
    return ResultCode::Success;
}

ResultCode aef::matrix::OneMklMatrixBackend::shutdown() {
    if (!ptr->init) {
        return ResultCode::S_NOTHING_PERFORMED;
    }
    delete ptr->device;
    delete ptr->queue;
    return ResultCode::Success;
}

ResultCode aef::matrix::OneMklMatrixBackend::set_max_size(int nMaxDim) {
    return ResultCode::Unimplemented;
}

ResultCode aef::matrix::OneMklMatrixBackend::multiply(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
    return ResultCode::Unimplemented;
}

ResultCode aef::matrix::OneMklMatrixBackend::commutator(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B, Eigen::MatrixXcd& out) {
    return ResultCode::Unimplemented;
}

ResultCode aef::matrix::OneMklMatrixBackend::group_action(Eigen::MatrixXcd& out, Eigen::MatrixXcd& U, Eigen::MatrixXcd& A) {
    return ResultCode::Unimplemented;
}

ResultCode aef::matrix::OneMklMatrixBackend::expectation_value(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A) {
    return ResultCode::Unimplemented;
}

ResultCode aef::matrix::OneMklMatrixBackend::matrix_element(dcomplex& out, Eigen::VectorXcd& v1, Eigen::MatrixXcd& A, Eigen::VectorXcd& v2) {
    return ResultCode::Unimplemented;
}

ResultCode aef::matrix::OneMklMatrixBackend::diagonalize(Eigen::MatrixXcd& mat, Eigen::VectorXcd& evals, Eigen::MatrixXcd& evecs) {
    return ResultCode::Unimplemented;
}

#pragma once
#include "aef/MatrixOpBackend.h"

namespace aef::matrix {
    struct OneMKLdata;
    class OneMklMatrixBackend : public IMatrixOpBackend {
        OneMKLdata* ptr;
    public:
        OneMklMatrixBackend();
        ~OneMklMatrixBackend();
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
};
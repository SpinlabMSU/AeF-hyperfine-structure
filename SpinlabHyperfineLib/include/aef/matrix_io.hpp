#pragma once
#include <string>
#include <Eigen/Eigen>
#include <iostream>
#include <fstream>
#include <zlib.h>
//#include <boost/iostreams/stream.hpp>


namespace Eigen {
    template<class Matrix>
    void write_binary(const char* filename, const Matrix& matrix) {
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
        out.write((char*)(&rows), sizeof(typename Matrix::Index));
        out.write((char*)(&cols), sizeof(typename Matrix::Index));
        out.write((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
        out.close();
    }
    template<class Matrix>
    void read_binary(const char* filename, Matrix& matrix) {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        typename Matrix::Index rows = 0, cols = 0;
        in.read((char*)(&rows), sizeof(typename Matrix::Index));
        in.read((char*)(&cols), sizeof(typename Matrix::Index));
        matrix.resize(rows, cols);
        in.read((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
        in.close();
    }

    template<typename... Args>
    void write_binary(std::ostream& out, const Eigen::DiagonalMatrix<Args...>& matrix) {
        typename Matrix::Index rows = matrix.rows(), cols = 1;
        out.write((char*)(&rows), sizeof(typename Matrix::Index));
        out.write((char*)(&cols), sizeof(typename Matrix::Index));
        out.write((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
    }

    template<class Matrix>
    void write_binary(std::ostream& out, const Matrix& matrix) {
        typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
        out.write((char*)(&rows), sizeof(typename Matrix::Index));
        out.write((char*)(&cols), sizeof(typename Matrix::Index));
        out.write((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
    }

    template<typename... Args>
    void read_binary(std::istream& in, Eigen::DiagonalMatrix<Args...>& matrix) {
        typename Matrix::Index rows = 0, cols = 1;
        in.read((char*)(&rows), sizeof(typename Matrix::Index));
        in.read((char*)(&cols), sizeof(typename Matrix::Index));
        matrix.resize(rows);
        in.read((char*)matrix.diagonal().data(), rows * cols * sizeof(typename Matrix::Scalar));
    }

    template<class Matrix>
    void read_binary(std::istream& in, Matrix& matrix) {
        typename Matrix::Index rows = 0, cols = 0;
        in.read((char*)(&rows), sizeof(typename Matrix::Index));
        in.read((char*)(&cols), sizeof(typename Matrix::Index));
        matrix.resize(rows, cols);
        std::cout << "Read " << rows << " " << cols << std::endl;
        in.read((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
    }



} // Eigen::

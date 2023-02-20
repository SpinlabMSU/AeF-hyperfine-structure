#pragma once
#include <string>
#include <Eigen/Eigen>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <format>
//#include <boost/iostreams/stream.hpp>
extern uint64_t stream_pos;

namespace Eigen {
#if 0
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

    template<typename... Args>
    void read_binary(std::istream& in, Eigen::DiagonalMatrix<Args...>& matrix) {
        typename Matrix::Index rows = 0, cols = 1;
        in.read((char*)(&rows), sizeof(typename Matrix::Index));
        in.read((char*)(&cols), sizeof(typename Matrix::Index));
        matrix.resize(rows);
        in.read((char*)matrix.diagonal().data(), rows * cols * sizeof(typename Matrix::Scalar));
    }

#endif
    
    template<class Matrix>
    void write_binary(std::ostream& out, const Matrix& matrix) {
        constexpr uint32_t magic = 0xffddeeff;
        out.write((char*)&magic, sizeof(magic));

        typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
        std::cout << std::format("Writing rows = {}, cols = {}", rows, cols) << std::endl;
        out.write((char*)(&rows), sizeof(typename Matrix::Index));
        out.write((char*)(&cols), sizeof(typename Matrix::Index));
        out.write((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));

        constexpr uint32_t emarker = 0xffddccff;
        out.write((char*)&emarker, sizeof(emarker));
    }


#ifndef AEF_STREAM_NDEBUG
    void print_stream_position(std::istream& in) {
        auto tellg = in.tellg();
        std::streamoff pos = tellg;
        std::cout << std::format("Stream position {} = 0x{:x}", pos, pos) << std::endl;
    }
#else
#define print_stream_position(in) ((void)in)
#endif
    template<class Matrix>
    void read_binary(std::istream& in, Matrix& matrix) {
        constexpr uint32_t magic = 0xffddeeff;
        print_stream_position(in);
        uint32_t rmagic = 0;
        in.read((char*)&rmagic, sizeof(rmagic));
        ::stream_pos += 4;
        print_stream_position(in);

        std::cout << std::format("Read magic {:x}", rmagic) << std::endl;

        if (magic != rmagic) {
#ifndef AEF_STREAM_NDEBUG
            for (int i = 1; i <= 32; i++) {
                in.read((char*)&rmagic, sizeof(rmagic));
                std::cout << std::format("Tried extra magic #{:x} : {:x}", i, rmagic) << std::endl;
            }
#endif
            DebugBreak();
            throw std::exception("BAD MAGIC");
        }

        typename Matrix::Index rows = 0, cols = 0;
        in.read((char*)(&rows), sizeof(typename Matrix::Index));
        ::stream_pos += sizeof(typename Matrix::Index);
        print_stream_position(in);
        
        in.read((char*)(&cols), sizeof(typename Matrix::Index));
        ::stream_pos += sizeof(typename Matrix::Index);
        print_stream_position(in);
        std::cout << "Read rows = " << rows << ", cols = " << cols << std::endl;
        matrix.resize(rows, cols);
        auto pre = in.tellg();

        size_t size = rows * cols * sizeof(typename Matrix::Scalar);

        in.read((char*)matrix.data(), size);
        ::stream_pos += rows * cols * sizeof(typename Matrix::Scalar);
#ifndef AEF_STREAM_NDEBUG
        print_stream_position(in);
        auto post = in.tellg();

        auto diff = post - pre;
        std::cout << std::format("Tried to read {} bytes, actually read {} bytes", size, diff) << std::endl;
#endif


        if (in.eof() || in.fail() || in.bad()) {
            std::cout << "Error eof=" << in.eof() << " fail=" << in.fail() << " bad=" << in.bad() << " tellg=" << in.tellg() << std::endl;
            //printf("ERROR eof=%d fail=%d bad=%d, %zu\n", in.eof(), in.fail(), in.bad(), in.tellg());
            throw std::exception("EOF or FAIL or BAD");
        }

        constexpr uint32_t emarker = 0xffddccff;
        uint32_t rmarker = 0;
        in.read((char*)&rmarker, sizeof(rmarker));
        ::stream_pos += sizeof(rmarker);
        print_stream_position(in);
        std::cout << std::format("Read end-marker {:x}", rmarker) << std::endl;

        if (rmarker != emarker) {
            throw std::exception("BAD END MARKER");
        }
    }



} // Eigen::

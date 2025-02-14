#pragma once

#include <Eigen/Eigen>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <zlib.h>
#include <fmt.hpp>


#include <TTree.h>

namespace aef {
    /// <summary>
    /// Write out a matrix representing a
    /// </summary>
    /// <param name="op">an Eigen::MatrixXcd representing a self-adjoint operator</param>
    /// <param name="name">The name of the TTree to save to</param>
    /// <param name="title">This will be passed as the title to the TTree constructor</param>
    /// <param name="mag_thresh"></param>
    /// <param name="use_rel_thresh"></param>
    /// <returns></returns>
    TTree* write_self_adj_op_tree(Eigen::MatrixXcd& op, const char* name, const char* title, double mag_thresh = 1.0e-17, bool use_rel_thresh = true);
    TTree* write_matrix_tree(Eigen::MatrixXcd& op, const char* name, const char* title, double mag_thresh = 1.0e-17, bool use_rel_thres = true, bool cartesian=true);
};
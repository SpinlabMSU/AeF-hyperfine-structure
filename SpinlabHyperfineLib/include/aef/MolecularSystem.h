#pragma once
#ifndef _AEF_MOLECULAR_SYSTEM_H

#define _AEF_MOLECULAR_SYSTEM_H 1
#include "Eigen/Eigen"
#include "aef.h"
#include "jf_basis_vec.h"
#include <filesystem>
#include <istream>
#include <unordered_map>
#include <vector>

namespace aef {
    /// <summary>
    /// Future V2 output
    /// </summary>
    class MolecularSystem {
        spin nmax;
        std::vector<jf_basis_vec> basis;
        std::unordered_map<std::string, Eigen::MatrixXcd*> operators;

    public:
        MolecularSystem();
        ~MolecularSystem();

        bool load(std::istream& in);
        bool save(std::ostream& out);

        // convienence methods
        bool load_matrix_elts(std::string inpath);
        bool save_matrix_elts(std::string out);

        bool load_matrix_elts(std::filesystem::path inpath);
        bool save_matrix_elts(std::filesystem::path out);
    };
};
#endif //_AEF_MOLECULAR_SYSTEM_H
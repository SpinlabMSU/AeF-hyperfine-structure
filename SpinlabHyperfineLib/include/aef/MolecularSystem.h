/*
    aef/MolecularSystem.h -- class that implements matrix-element calculations
    for the 225RaF molecule.

    This file is part of the AeF-hyperfine-structure program. 
    
    AeF-hyperfine-structure is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version.

    AeF-hyperfine-structure is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    AeF-hyperfine-structure. If not, see <https://www.gnu.org/licenses/>.
*/
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

        struct parameters {
            
        };

    public:
        MolecularSystem(spin nmax_ = 0.0, double E_z = 1.0, double K = 0.0);
        ~MolecularSystem();

        void set_nmax(spin nmax_);


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
/*
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
#include "pch.h"
#include "aef/MolecularSystem.h"

namespace aef {
    MolecularSystem::MolecularSystem(spin nmax_, double E_z, double K):
        nmax(nmax_)
    {
    }
    void MolecularSystem::set_nmax(spin nmax_) {
    }
    MolecularSystem::~MolecularSystem() {
    }
    bool MolecularSystem::load(std::istream& in) {
        return false;
    }
    bool MolecularSystem::save(std::ostream& out) {
        return false;
    }
};
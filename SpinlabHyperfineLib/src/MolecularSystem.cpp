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
#include <unordered_map>

namespace aef {
    MolecularSystem::MolecularSystem(IMolecularCalculator &calc_, spin nmax_, double E_z_, double K_) :
        calc(calc_), E_z(E_z_), K(K_) {
        diagonalized = false;
        init = false;
        dkq_init = false;
        enableDev = K != 0;

        set_nmax(nmax_);
    }
    void MolecularSystem::set_nmax(spin nmax_) {
        // nmax must be a non-negative half-integer
        assert(nmax >= 0 && floor(nmax * 2) == (nmax * 2));
        nmax = nmax_;
        calc.set_nmax(nmax);
        nBasisElts = calc.get_nBasisElts();
        H_rot.resize(nBasisElts);
        H_hfs.resize(nBasisElts, nBasisElts);
        H_stk.resize(nBasisElts, nBasisElts);
        H_dev.resize(nBasisElts, nBasisElts);
        H_tot.resize(nBasisElts, nBasisElts);
        F_z.resize(nBasisElts, nBasisElts);
        init = false;
        d1t.resize(nBasisElts, nBasisElts);
        d10.resize(nBasisElts, nBasisElts);
        d11.resize(nBasisElts, nBasisElts);
        dkq_init = false;
        Es.resize(nBasisElts);
        Vs.resize(nBasisElts, nBasisElts);
        Vst.resize(nBasisElts, nBasisElts);
        diagonalized = false;

        for (auto& [id, mat] : opMatMap) {
            mat->resize(nBasisElts, nBasisElts);
        }

    }
    void MolecularSystem::calculate_matrix_elts() {
        calc.calculate_H_rot(this->H_rot);
        calc.calculate_H_hfs(this->H_hfs);
        calc.calculate_H_stk(this->H_stk);
        calc.calculate_H_dev(this->H_dev);

        H_tot.setZero();
        H_tot.diagonal() = H_rot.diagonal();
        H_tot += H_hfs + E_z * H_hfs + K * H_dev;

        init = true;
    }
    void MolecularSystem::calculate_dkq() {
        calc.calculate_dkq(this->d1t, -1);
        calc.calculate_dkq(this->d10,  0);
        calc.calculate_dkq(this->d11, +1);
        dkq_init = true;
    }



    MolecularSystem::~MolecularSystem() {
        this->H_rot.resize(0);
        this->H_hfs.resize(0, 0);
        this->H_stk.resize(0, 0);
        this->H_dev.resize(0, 0);
        this->H_tot.resize(0, 0);
        this->d1t.resize(0, 0);
        this->d10.resize(0, 0);
        this->d11.resize(0, 0);

        for (auto& [id, mat] : opMatMap) {
            delete mat;
        }
        opMatMap.clear();
        opMap.clear();
    }

    aef::ResultCode MolecularSystem::diagonalize() {
        //using aef::ResultCode;
        auto rc = aef::ResultCode::Success;
        rc = aef::matrix::set_max_size(nBasisElts);
        if (aef::failed(rc)) return rc;
        if (enableDev) {
            rc = aef::matrix::diagonalize(H_tot, Es, Vs);
        } else {
            // Free-space: m_f is a good quantum number, want to simultaneously
            // diagonalize H_tot and F_z This is done using the method given in
            // https://math.stackexchange.com/a/4388322 and proven in
            // https://math.stackexchange.com/a/3951339.  However, I'm omitting the
            // randomization portion because it shouldn't be necessary (if there's some
            // small mixing of the m_f it doesn't really matter, and in practice they
            // don't mix.
            constexpr dcomplex t = 100.0; // +15i;
            // naughty hack: H_dev doesn't actually contain anything when enableDev == false, so we can use it
            // as our temporary here instead of making a new temporary matrix
            H_dev = H_tot + t * F_z;
            rc = aef::matrix::diagonalize(H_dev, Es, Vs);
            rc = aef::matrix::group_action(H_dev, Vs, H_tot);
            Es = H_dev.diagonal();
            H_dev.setZero();
        }
        diagonalized = true;
        return rc;
    }

    /////// IO code
     // load and save --> todo choose either XIFF or TTree
    aef::ResultCode MolecularSystem::load(std::istream& in) {
        // XIFF draft -->
        // XIFF header type "AEF0"
        

        return aef::ResultCode::Unimplemented;
    }
    aef::ResultCode MolecularSystem::save(std::ostream& out) {
        return aef::ResultCode::Unimplemented;
    }

    aef::ResultCode MolecularSystem::load(std::string inpath) {
        std::ifstream in(inpath, std::ios::binary);
        return load(in);
    }

    aef::ResultCode MolecularSystem::save(std::string outpath) {
        std::ofstream out(outpath, std::ios::binary | std::ios::trunc | std::ios::out);
        return save(outpath);
    }

    aef::ResultCode MolecularSystem::load(std::filesystem::path inpath) {
        std::ifstream in(inpath, std::ios::binary);
        return load(in);
    }

    aef::ResultCode MolecularSystem::save(std::filesystem::path outpath) {
        std::ofstream out(outpath, std::ios::binary | std::ios::trunc | std::ios::out);
        return save(outpath);
    }


    ////// PTFW2

    aef::operators::IOperator* MolecularSystem::getOperator(const std::string& id) {
        if (!opMap.contains(id)) {
            return nullptr;
        }
        return opMap[id];
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="id"></param>
    /// <param name="op"></param>
    void MolecularSystem::addOperator(const std::string& id, aef::operators::IOperator* op) {
        opMap[id] = op;
        opMatMap[id] = new Eigen::MatrixXcd(nBasisElts, nBasisElts); // allocate matrix
    }

    void MolecularSystem::evaluate(void) {
        for (auto& [id, op] : opMap) {
            Eigen::MatrixXcd* opMat;
            if (!opMatMap.contains(id)) {
                opMat = new Eigen::MatrixXcd(nBasisElts, nBasisElts);
                opMatMap.emplace(id, opMat);
            } else {
                opMat = opMatMap[id];
            }
            op->fillMatrix(*opMat);
        }
    }

    Eigen::MatrixXcd* MolecularSystem::getOperatorMatrix(const std::string& id) {
        return opMatMap[id];
    }

    dcomplex MolecularSystem::get_matrix_element(const std::string& id, int eidx1, int eidx2) {
        Eigen::VectorXcd state1 = Vs.col(eidx1);
        Eigen::VectorXcd state2 = Vs.col(eidx2);
        Eigen::MatrixXcd* op = getOperatorMatrix(id);

        dcomplex out;
        aef::matrix::matrix_element(out, state1, *op, state2);
        return out;
    }

    dcomplex MolecularSystem::expectation_value(const std::string& id, int eidx1) {
        Eigen::VectorXcd state1 = Vs.col(eidx1);
        Eigen::MatrixXcd* op = getOperatorMatrix(id);

        dcomplex out;
        aef::matrix::expectation_value(out, state1, *op);
        return out;
    }

    aef::ResultCode MolecularSystem::delta_E_lo(const std::string& id, Eigen::VectorXcd& output, Eigen::MatrixXcd* workspace) {
        bool internal_workspace = !workspace;
        if (internal_workspace) {
            workspace = new Eigen::MatrixXcd;
            workspace->resize(nBasisElts, nBasisElts);
        }
        Eigen::MatrixXcd* op = getOperatorMatrix(id);
        // get expectation values
        aef::matrix::group_action(*workspace, Vs, *op);
        output = workspace->diagonal();
        // second order ???

        if (internal_workspace) {
            delete workspace;
        }
    }

    /// Makermap code
    namespace {
        std::unordered_map<std::string, IMolecularCalculator::pfnMolCalcMaker> makerMap;
    };

    void IMolecularCalculator::registerMolCalcType(std::string name, pfnMolCalcMaker ctor) {
        makerMap[name] =  ctor;
    }
    IMolecularCalculator* IMolecularCalculator::makeCalculatorOfType(std::string name) {
        if (makerMap.contains(name)) {
            auto ctor = makerMap[name];
            return ctor();
        }
        return nullptr;
    }
};
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
#include <aef/xiffstream.h>
#include <aef/io/aefchunk.h>

#include <TKey.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TMap.h>

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

    struct prms_payload_fixed {
        uint32_t twice_nmax;
        int n_extra;
        double E_z;
        double K;
    };

    namespace molsys_flags {
        constexpr uint16_t compress = 0x1;
    };

    namespace {
        /*
         * <summary>
         * This enum class describes what each bit of flags means.
         * </summary>
         */
        enum hyp_flags : uint16_t {
            // 
            FLAG_INIT = 1,
            // if this flag is set, then the file contains the computed eigenenergies/eigenstates of the Hamiltonian
            FLAG_DIAG = 2,
            // if this flag is set, then the file contains the MDA spherical tensor operators
            FLAG_DKQ_INIT = 4,
            // if this flag is set, then the file contains the specific values of electric field strength and 
            // devonshire coupling constant used to calculate these parts of the Hamiltonian.
            FLAG_E_DEV_PARMS = 8
        };
    };

    aef::ResultCode MolecularSystem::read_chunk(std::istream& in, void* chdr_, void *dst) {
        auto *pChdr = static_cast<aef::chunk::chunk_hdr*>(chdr_);
        auto &chdr = *pChdr;
        in.read((char*)&chdr, sizeof(chdr));
        switch (chdr.type) {
        case aef::chunk::hamiltonian_chunk:
            // hamiltonian_chunk
            std::cout << "[aef::MolecularSystem] reading hamiltonian chunk" << std::endl;
            aef::chunk::chunk_hdr mhdr;
            read_chunk(in, (void*)&mhdr, (void*)&H_rot.diagonal());
            read_chunk(in, (void*)&mhdr, (void*)&H_hfs);
            read_chunk(in, (void*)&mhdr, (void*)&H_stk);
            read_chunk(in, (void*)&mhdr, (void*)&H_dev);
            break;
        case aef::chunk::matrix:
            if (dst == nullptr) {
                std::clog << "[aef::MolecularSystem] read_chunk error, no destination provided for matrix chunk" << std::endl;
            }
            break;
        default:
            std::clog << "[aef::MolecularSystem] Warning: encountered unknown tag during loading" << std::endl;
            // handle end tag
        case aef::chunk::end0:
            std::cout << "[aef::MolecularSystem] load complete, hit end tag" << std::endl;
            break;
        }

        return aef::ResultCode();
    }


    /////// IO code
     // load and save --> todo choose either XIFF or TTree
    aef::ResultCode MolecularSystem::load(std::istream& in_, char *path) {
        std::cout << fmt::format("[aef::MolecularSystem] loading from stream, path {}", path) << std::endl;
        ////////////////
        //// aefchunk file format
        
        /// file header
        aef::chunk::file_hdr fhdr = {};
        in_.read((char*)&fhdr, sizeof(fhdr));
        
        // test magic
        if (fhdr.hdr.type != aef::chunk::file_magic) {
            std::clog << fmt::format("[MolecularSystem aefchunk loader] Error: file \"{}\" is not an aefchunk file!", path) << std::endl;
            return aef::ResultCode::InvalidFormat;
        }

        // test filetype
        if (fhdr.filetype != aef::chunk::mol_sys) {
            std::string str((char*)&(fhdr.filetype), 4);
            std::clog << fmt::format("[MolecularSystem aefchunk loader] Error: file \"{}\" has invalid format \"{}\"!", path, str) << std::endl;
            return aef::ResultCode::InvalidFormat;
        }
        // check version is in-range
        uint16_t version = fhdr.hdr.version;
        if (version < (uint16_t)MINIMUM_LOAD_VERSION) {
            std::cout << "" << std::endl;
        }
        
        /// flags
        uint16_t flags = fhdr.hdr.flags;
        // handle compression.
        zstr::istream zin(in_);
        bool is_stream_compressed = (flags & (uint16_t)molsys_flags::compress);
        std::istream* pIn = is_stream_compressed ? new zstr::istream(in_) : &in_;
        std::istream& in = *pIn;

        /// Start reading chunks
        aef::chunk::chunk_hdr chdr;

        // the params block must come first
        in.read((char*)&chdr, sizeof(chdr));
        if (chdr.type != aef::chunk::prms) {
            std::clog << "[aef::MolecularSystem] Error: molsys file {} is malformed" << std::endl;
        }
        this->init = flags & FLAG_INIT;
        this->dkq_init = flags & FLAG_DKQ_INIT;
        this->diagonalized = flags & FLAG_DIAG;

        {
            // handle params chunk "payload".  Here the "version" is actually used to store the payload size
            uint16_t payload_size = chdr.version;
            prms_payload_fixed *pay = (prms_payload_fixed*)calloc(payload_size, 1);
            in.read((char*)pay, payload_size);
            this->nmax = pay->twice_nmax / 2.0;

            this->set_nmax(nmax);
            this->E_z = pay->E_z;
            this->K = pay->K; 
            // todo add


            free(pay);
        }
        // handle the other chunks in-order using read_chunk
        do {
            this->read_chunk(in, (void*)&chdr, nullptr);
        } while (chdr.type != aef::chunk::end0);

        // finish
        if (pIn != &in_) {
            delete pIn; pIn = nullptr;
        }
        return aef::ResultCode::Success;
        /////////////////////////////////////
        //// XIFF code -- format aef0
        //// note that file header will be uncompressed
        xiff::xiff_hdr xhdr = {};
        in.read((char*) & xhdr, sizeof(xhdr));
        
        // match type
        if (xhdr.file_hdr.type != xiff::common_cc::xiff || xhdr.file_hdr.version > 0) {
            std::clog << std::format("[aef::MolecularSystem] file {} is not a recognized XIFF file", path) << std::endl;
            return ResultCode::InvalidFormat;
        }

        // match file type
        if (xhdr.ftype != xiff::aef_cc::aef0) {
            std::clog << std::format("[aef::MolecularSystem] XIFF file {}  is not an AeFDat file", path) << std::endl;
            return ResultCode::InvalidFormat;
        }

        xiff::chunk_hdr chdr = {};
        in.read((char*)&chdr, sizeof(chdr));
        if (chdr.type != xiff::aef_cc::prms) {
            return ResultCode::InvalidFormat;
        }

        do {
            chdr = {};
            in.read((char*)&chdr, sizeof(chdr));
            switch (std::bit_cast<uint32_t>(chdr.type.cc)) {
            case std::bit_cast<uint32_t>(xiff::aef_cc::oplist.cc):
                break;
            }
        } while (chdr.type != xiff::common_cc::end0);


        return aef::ResultCode::Unimplemented;
    }
    aef::ResultCode MolecularSystem::save(std::ostream& out_, char *path) {
        return aef::ResultCode::Unimplemented;
    }

    aef::ResultCode MolecularSystem::load(std::string inpath) {
        TFile f(inpath.c_str(), "READ");

        std::ifstream in(inpath, std::ios::binary);
        return load(in);
    }

    enum class aef_molsys_flags:uint64_t {
        flag_init = 1,
        flag_diag = 2,
        flag_dkq_init = 4,
        flag_enableDev = 8,
    };

    aef::ResultCode MolecularSystem::save(std::string outpath) {
        TFile sav(outpath.c_str(), "RECREATE");
        sav.cd();
        // save properties: nmax, E_z, K, 
        TParameter<uint64_t> prop_version("version", 0);
        uint64_t flags = 0;
        if (this->init) flags |= (uint64_t)aef_molsys_flags::flag_init;
        if (this->diagonalized) flags |= (uint64_t)aef_molsys_flags::flag_diag;
        if (this->dkq_init) flags |= (uint64_t)aef_molsys_flags::flag_dkq_init;
        if (this->enableDev) flags |= (uint64_t)aef_molsys_flags::flag_enableDev;

        // bitset
        TParameter<uint64_t> prop_flags("flags", flags);
        TParameter<spin> prop_nmax("nmax", nmax);
        TParameter<double> prop_E_z("E_z", E_z);
        TParameter<double> prop_K("K", K);


        TMap m;
        m.Add(&prop_version);
        m.Add(&prop_nmax);
        m.Add(&prop_E_z);
        m.Add(&prop_K);
        m.Write("properties", TObject::kSingleKey);
        // save operator list

        //


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
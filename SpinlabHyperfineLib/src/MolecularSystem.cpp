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
#include <aef/matrix_io.hpp>

#include <TKey.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TMap.h>
#include <spanstream>

namespace aef {
    MolecularSystem::MolecularSystem(IMolecularCalculator *calc_, spin nmax_, double E_z_, double K_) :
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
        calc->set_nmax(nmax);
        nBasisElts = calc->get_nBasisElts();
        std::cout << fmt::format("Resizing aef::MolecularSystem to nmax={}, will have {} basis elements.", nmax_, nBasisElts) << std::endl;
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
        calc->calculate_H_rot(this->H_rot);
        calc->calculate_H_hfs(this->H_hfs);
        calc->calculate_H_stk(this->H_stk, E_z);
        calc->calculate_H_dev(this->H_dev, K);

        H_tot.setZero();
        H_tot.diagonal() = H_rot.diagonal();
        H_tot += H_hfs + H_stk + H_dev;

        init = true;

        calc->calculate_dkq(this->d1t, -1);
        calc->calculate_dkq(this->d10, 0);
        calc->calculate_dkq(this->d11, +1);
        dkq_init = true;
    }

    aef::IMolecularCalculator* MolecularSystem::get_calc() {
        return calc;
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

        uint32_t cbCalcData;
        uint32_t calcTypeLen;
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
        };
    };

    aef::ResultCode MolecularSystem::write_chunk(std::ostream& out, void* chdr, void* data) {
        return aef::ResultCode();
    }

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
        {
            if (dst == nullptr) {
                std::clog << "[aef::MolecularSystem] read_chunk error, no destination provided for matrix chunk" << std::endl;
                return aef::ResultCode::InvalidArgument;
            }

            ptrdiff_t offset = chdr.version | (chdr.flags & 0xff) << 8;
            char* dst_ = ((char*)this) + offset;
            if (chdr.flags & (1 << 15)) {
                dst = dst_;
            }
            // flags indicates the datatype
            bool is_vector = chdr.flags & 1;
            bool is_complex = chdr.flags & 2;

            if (is_complex) {
                if (is_vector) {
                    auto* v = (Eigen::VectorXcd*)dst;
                    Eigen::read_binary(in, *v);
                } else {
                    auto* v = (Eigen::MatrixXcd*)dst;
                    Eigen::read_binary(in, *v); 
                }
            }
            return aef::ResultCode::Success;
        }
        case aef::chunk::oplist:
            // note: oplist stores version in upper byte of flags
            std::clog << "[aef::MolecularSystem] loading perturbative operator map" << std::endl;
            {
                std::vector<std::string> ops;
                int version = chdr.flags >> 8;
                int nops = (chdr.flags & 0xff) | chdr.version;
                char opId[max_op_id_len];
                memset(opId, 0, max_op_id_len);
                for (int idxOp = 0; idxOp < nops; idxOp++) {
                    int len = 0;
                    in.read((char*)&len, sizeof(len));
                    assert(len < max_op_id_len);
                    in.read(opId, len);
                    // TODO figure out how to load IOperator
                }
                //for (int idx)
            };
            break;
        default:
            std::clog << "[aef::MolecularSystem] Warning: encountered unknown tag during loading" << std::endl;
            // handle end tag
            [[fallthrough]];
        case aef::chunk::end0:
            std::cout << "[aef::MolecularSystem] load complete, hit end tag" << std::endl;
            return aef::ResultCode::Success;
        }

        return aef::ResultCode();
    }


    /////// IO code
    // aefchunk load code
    aef::ResultCode MolecularSystem::load(std::istream& in_, const char *path) {
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
            std::clog << fmt::format("[MolecularSystem aefchunk loader] Error file \"{}\" uses unsupported version {} (too old)", path, version) << std::endl;
            return aef::ResultCode::InvalidFormat;
        }

        if (version > (uint16_t)MAXIMUM_LOAD_VERSION) {
            std::clog << fmt::format("[MolecularSystem aefchunk loader] Error file \"{}\" uses unsupported version {} (too old)", path, version) << std::endl;
            return aef::ResultCode::InvalidFormat;
        }
        
        /// flags
        uint16_t flags = fhdr.hdr.flags;
        // handle compression.
        zstr::istream zin(in_);
        bool is_stream_compressed = (flags & (uint16_t)molsys_flags::compress);
        std::istream* pIn = is_stream_compressed ? new zstr::istream(in_) : &in_;
        std::istream& in = *pIn;

        /// Start reading chunks
        aef::chunk::chunk_hdr chdr = {};

        // the params block must come first
        in.read((char*)&chdr, sizeof(chdr));
        if (chdr.type != aef::chunk::prms) {
            std::clog << "[aef::MolecularSystem] Error: molsys file {} is malformed" << std::endl;
            return aef::ResultCode::InvalidFormat;
        }
        this->init = flags & FLAG_INIT;
        this->dkq_init = flags & FLAG_DKQ_INIT;
        this->diagonalized = flags & FLAG_DIAG;

        {
            // handle params chunk "payload".  Here the "version" is actually used to store the payload size
            uint16_t payload_size = chdr.version;
            prms_payload_fixed *pay = (prms_payload_fixed*)calloc(payload_size, 1);
            assert(pay);
            in.read((char*)pay, payload_size);
            this->nmax = pay->twice_nmax / 2.0;

            this->set_nmax(nmax);
            this->E_z = pay->E_z;
            this->K = pay->K;

            //// read aef::IMolecularCalculator
            char* curr = ((char*)pay) + sizeof(prms_payload_fixed);
            std::string calcType(curr, pay->calcTypeLen);
            curr += pay->calcTypeLen;

            // use a seperate input stream to try to prevent 
            std::ispanstream calcIn(std::span(curr, pay->cbCalcData));

            this->calc = aef::IMolecularCalculator::makeCalculatorOfType(calcType);
            if (calc) {
                this->calc->set_nmax(nmax);
                calc->load(calcIn);
            }
            free(pay);
            if (!calc) {
                std::clog << fmt::format("[aef::MolecularSystem] Error: molsys file \"{}\" uses unavailable IMolecularCalculator type \"{}\".",
                    path, calcType)<< std::endl;
                return aef::ResultCode::NotAvailable;
            }
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
    }



    aef::ResultCode MolecularSystem::save(std::ostream& out_, const char *path) {
        // TODO implement
        // write file header
        std::ostream* out = &out_;
        aef::chunk::file_hdr hdr = {};
        {
            // always compress 
            uint16_t flags = (uint16_t)molsys_flags::compress;
            hdr.hdr = { .type = aef::chunk::file_magic, .version = (uint16_t)CURRENT_SAVE_VERSION, .flags = flags };
            hdr.filetype = aef::chunk::mol_sys;
        }
        out_.write((char*)&hdr, sizeof(hdr));

        if (hdr.hdr.flags & (uint16_t)molsys_flags::compress) {
            out = new zstr::ostream(out_);
        }

        // write parameter chunk
        aef::chunk::chunk_hdr chdr = {};
        {
            chdr.type = aef::chunk::prms;
            chdr.version = 0;
            chdr.flags = 0;
            if (this->init) chdr.flags |= FLAG_INIT;
            if (this->diagonalized) chdr.flags |= FLAG_DIAG;
            if (this->dkq_init) chdr.flags |= FLAG_DKQ_INIT;
            std::string calcType = calc->get_calc_type();
            struct prms_payload_fixed pay = {};
            pay.calcTypeLen = (int)calcType.length();

            std::ostringstream os;
            calc->save(os);
            auto buf = os.str();

            out->write((char*)&chdr, sizeof(chdr));
            out->write((char*)&pay, sizeof(pay));
            out->write(calcType.data(), pay.calcTypeLen);
            out->write(buf.data(), buf.length());
        }
        // write matricies
        write_vector(*out, &this->H_rot.diagonal());
        write_matrix(*out, &this->H_hfs);
        // write operators
        return aef::ResultCode::Unimplemented;
    }

    aef::ResultCode aef::MolecularSystem::write_matrix(std::ostream& out, Eigen::MatrixXcd* mat) {
        return aef::ResultCode::Success;
    }

    aef::ResultCode aef::MolecularSystem::write_vector(std::ostream& out, Eigen::VectorXcd* vec) {
        return aef::ResultCode::Success;
    }

    aef::ResultCode MolecularSystem::load(std::string inpath) {
        std::ifstream in(inpath, std::ios::binary);
        return load(in, inpath.c_str());
    }

    enum class aef_molsys_flags:uint64_t {
        flag_init = 1,
        flag_diag = 2,
        flag_dkq_init = 4,
        flag_enableDev = 8,
    };

    aef::ResultCode MolecularSystem::save(std::string outpath) {
        std::ofstream out(outpath, std::ios::binary | std::ios::trunc | std::ios::out);
        return save(out, outpath.c_str());
    }

    aef::ResultCode MolecularSystem::load(std::filesystem::path inpath) {
        std::ifstream in(inpath, std::ios::binary);
        return load(in);
    }

    aef::ResultCode MolecularSystem::save(std::filesystem::path outpath) {
        std::ofstream out(outpath, std::ios::binary | std::ios::trunc | std::ios::out);
        return save(out, outpath.string().c_str());
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
        aef::ResultCode rc = ResultCode::S_NOTHING_PERFORMED;

        if (internal_workspace) {
            workspace = new Eigen::MatrixXcd;
            workspace->resize(nBasisElts, nBasisElts);
        }
        Eigen::MatrixXcd* op = getOperatorMatrix(id);
        
        if (!op) {
            rc = ResultCode::InvalidArgument;
            goto end;
        }

        // get expectation values
        rc = aef::matrix::group_action(*workspace, Vs, *op);
        output = workspace->diagonal();
        rc = aef::ResultCode::Success;
        // second order ??? 

        end:
        if (internal_workspace) {
            delete workspace;
        }
        return rc;
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
    std::string universal_diatomic_basis_vec::ket_string() {
        if (type == coupling_type::j_basis) {
            return fmt::format("|n={}, j={}, f={}, m_f={}>", n, j, f, m_f);
        } else if (type == coupling_type::jf_basis) {
            return fmt::format("|n={}, j={}, f1={}, f={}, m_f={}>", n, j, f_1, f, m_f);
        }
        return fmt::format("|i1={}, i2={}, s={}, l={}, r={}, n={}, j={}, f1={}, f={}, m_f={}>", 
            i1, i2, s, l, r,
            n, j, f_1, f, m_f);
    }
    std::string universal_diatomic_basis_vec::ket_csv_str() {
        if (type == coupling_type::j_basis) {
            return fmt::format("|n={} j={} f={} m_f={}>", n, j, f, m_f);
        } else if (type == coupling_type::jf_basis) {
            return fmt::format("|n={} j={} f1={} f={} m_f={}>", n, j, f_1, f, m_f);
        }
        return fmt::format("|i1={} i2={} s={} l={} r={} n={} j={} f1={} f={} m_f={}>",
            i1, i2, s, l, r,
            n, j, f_1, f, m_f);
    }
};

std::ostream& operator<<(std::ostream& os, universal_diatomic_basis_vec& v) {
    return (os << fmt::format("{}", v));
}

bool operator==(const universal_diatomic_basis_vec& v1, const universal_diatomic_basis_vec& v2) {
    if (v1.type != v2.type) return false;
    if (v1.i1 != v2.i1) return false;
    if (v1.i2 != v2.i2) return false;
    if (v1.s != v2.s) return false;
    if (v1.l != v2.l) return false;
    if (v1.r != v2.r) return false;
    if (v1.n != v2.n) return false;
    if (v1.j != v2.j) return false;
    if (v1.f_1 != v2.f_1) return false;
    if (v1.f != v2.f) return false;
    if (v1.m_f != v2.m_f) return false;

    return true;
}

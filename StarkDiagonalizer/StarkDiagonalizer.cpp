// StarkDiagonalizer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/*
    StarkDiagonalizer/StarkDiagonalizer.cpp -- this program diagonalizes the stark hamiltonian
    using the same j-basis that the main AeF-hyperfine-structure program uses.

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
#include <pch.h>
#include <aef/aef.h>
#include <aef/debug_stream.h>
#include <aef/matrix_utils.h>
#include <aef/teestream.hpp>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <cxxopts.hpp>
#include <fmt.hpp>
#include <aef/quantum.h>

namespace fs = std::filesystem;
using namespace std::chrono;
using aef::log_time_at_point;
using namespace aef::quantum;

#include "../AeF-hyperfine-structure.inl"


// the number of states in a singlet-triplet group pair (f=0 and f=1)
constexpr int num_singlet_triplet = 4;
// +-X +-Y +-Z
constexpr int num_orientations = 6;
// the total number of states in the lowest lying (in energy) "bottom group" of states
constexpr int bottom_group_size = num_singlet_triplet * num_orientations;

// starting/base index of +Z oriented states
constexpr int bidx_posz = 0;
// starting/base index of -Z oriented states
constexpr int bidx_negz = bottom_group_size - num_singlet_triplet;
// start of +- XY oriented states -- don't know how to distinguish which is which at the present time
constexpr int bidx_pmxy = num_singlet_triplet;

void dump_state(HyperfineCalculator& calc, int state, std::ostream& out) {
    out << state;
    // BUGFIX: states are column vectors not row vectors
    Eigen::VectorXcd state_vec = calc.Vs.col(state);
    for (int kdx = 0; kdx < calc.nBasisElts; kdx++) {
        dcomplex ampl = state_vec[kdx];
        out << "," << std::real(ampl) << "," << std::imag(ampl);
    }
    out << std::endl;
}

/// <summary>
/// Outputs some key information about each energy eigenstate including:
/// * the MDA expectation value (re then IM)
/// * the expectation values of n, j, f, and m_f
/// </summary>
/// <param name="output">the stream to output to</param>
/// <param name="calc"></param>
void output_state_info(std::ostream& output, HyperfineCalculator& calc
#ifndef DONT_USE_CUDA
    , Eigen::MatrixXcd& vals
#endif
) {
    output << "Index n, Energy (MHz), Re(<n|dx|n>), Re(<n|dy|n>), Re(<n|dz|n>), "
        "Im(<n|dx|n>), Im(<n|dy|n>), Im(<n|dz|n>), "
        "Re(<n|n|n>), Re(<n|j|n>), Re(<n|f|n>), Re(<n|m_f|n>),"
        "Im(<n|n|n>), Im(<n|j|n>), Im(<n|f|n>), Im(<n|m_f|n>),"
        "<n|(-1)^n|n>"
        << std::endl;

#ifndef DONT_USE_CUDA
    Eigen::VectorXcd d10s;
    Eigen::VectorXcd d11s;
    Eigen::VectorXcd d1ts;
    aef::matrix::group_action(vals, calc.Vs, calc.d10);
    d10s = vals.diagonal();
    aef::matrix::group_action(vals, calc.Vs, calc.d11);
    d11s = vals.diagonal();
    aef::matrix::group_action(vals, calc.Vs, calc.d1t);
    d1ts = vals.diagonal();
#endif
    for (size_t n = 0; n < calc.nBasisElts; n++) {
#ifdef DONT_USE_CUDA
        auto e_n = calc.Vs.col(n);
        // molecular dipole vector in spherical tensor form
        dcomplex d10 = expectation_value(e_n, calc.d10);
        dcomplex d11 = expectation_value(e_n, calc.d11);
        dcomplex d1t = expectation_value(e_n, calc.d1t);
#else
        dcomplex d10 = d10s(n);
        dcomplex d11 = d11s(n);
        dcomplex d1t = d1ts(n);
#endif
        constexpr double inv_sqrt2 = std::numbers::sqrt2 / 2.0;

        // convert to cartesian
        using namespace std::complex_literals;
        dcomplex dx = (d1t - d11) * inv_sqrt2;
        dcomplex dy = (d1t + d11) * 1i * inv_sqrt2;
        dcomplex dz = d10;

        // basis operator expectation values
        j_basis_vec v = expectation_values_jsq(calc, n);

        // output
        auto mda_ifo =
            fmt::format("{}, {}, {}, {}, {}, {}, {}, {}", n, std::real(calc.Es[n]),
                std::real(dx), std::real(dy), std::real(dz), std::real(dx),
                std::imag(dy), std::imag(dz));
        auto re_njfmf = fmt::format("{},{},{},{}", std::real(v.n), std::real(v.j), std::real(v.f), std::real(v.m_f));
        auto im_njfmf = fmt::format("{},{},{},{}", std::imag(v.n), std::imag(v.j), std::imag(v.f), std::imag(v.m_f));
        output << mda_ifo << ", " << re_njfmf << ", " << im_njfmf << "," << expect_parity(calc, n) << std::endl;
    }
    output.flush();
}

int main(int argc, char **argv) {
    std::cout << "Hello World!\n";
    auto dpath = fs::path("output");
    std::chrono::time_point<std::chrono::system_clock> start_time =
        std::chrono::system_clock::now();
    std::string stime = fmt::format("{0:%F}-{0:%H%M}{0:%S}", start_time);
    std::chrono::time_point<std::chrono::system_clock> prev_time = start_time;
    dpath /= stime;
    fs::create_directories(dpath);

    /// <summary>
    /// This is the field strength (in V/cm) used to calculate H_stk.
    /// DO NOT CHANGE THIS.
    /// 50 kV/cm = 25170 MHz/D is a standardized value used across the AeF-hyperfine-structure codebase
    /// </summary>
    constexpr double calc_E_z_V_per_cm = 50 * 1000;

    /// <summary>
    /// 50 kV/cm = 25170 MHz/D is the field strength used to calculate H_stk.
    /// DO NOT CHANGE THIS.
    /// This is a standardized value used across the AeF-hyperfine-structure codebase
    /// </summary>
    // note: can't substitute * 50 * 1000 with calc_E_z_V_per_cm because it rounds differently, giving
    // 25170.000000000004 isntead of 25170.0, which messes up file names
    constexpr double calc_E_z = unit_conversion::MHz_D_per_V_cm * 50 * 1000;

    // parameter block
    int param_nmax = 20;
    bool enable_debug_log = false;
    bool load_from_file = false;
    std::string loadname = "";
    bool print_extras = true;
    bool output_Es = true;

    // arg parsing
    cxxopts::Options options("stark_diagonalizer", "Program to simulate the hyperfine structure of"
        "diatomic Alkaline - monoflouride molecules");
    options.add_options()
        ("h,help", "Print usage")
        ("n,n_max", "Maximum n level to include", cxxopts::value<int>())
        ("d,enable_debug", "Enable debug mode", cxxopts::value<bool>()->default_value("false"))
        ("print_extras", "Print extra information", cxxopts::value<bool>()->default_value("true"))
        ("l,load", "Load molecular system operators from file", cxxopts::value<std::string>());

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("enable_debug")) {
        enable_debug_log = result["enable_debug"].as<bool>();
    }

    if (result.count("n_max")) {
        param_nmax = result["n_max"].as<int>();
    }
    if (result.count("load")) {
        load_from_file = true;
        loadname = result["load"].as<std::string>();
    }

    if (result.count("print_extras")) {
        print_extras = result["print_extras"].as<bool>();
    }

    fs::path odir = dpath;
    {
        std::ofstream ftyp(odir / "STARK_ONLY.txt");
    }

    // create info log
    std::ofstream oLog;

    debug_stream::debug_ostream* pDebug = new debug_stream::debug_ostream;
    teestream::teebuf* pBuf, * pOutb, * pErrb;

    std::streambuf* pLogBuf;
    std::streambuf* orig_coutb = std::cout.rdbuf();
    std::streambuf* orig_cerrb = std::cerr.rdbuf();

    {
        auto fpath = dpath / "out.log";
        oLog = std::ofstream(fpath, std::ios::trunc | std::ios::out);
        // test if logging to debug window enabled
        if (enable_debug_log) {
            // if enabled, tee to both the log file and the debug window
            pBuf = new teestream::teebuf(oLog.rdbuf(), pDebug->rdbuf());
            pLogBuf = pBuf;
        } else {
            // otherwise just output to the log file
            pBuf = nullptr;
            pLogBuf = oLog.rdbuf();
        }
        pOutb = new teestream::teebuf(pLogBuf, orig_coutb);
        std::cout.rdbuf(pOutb);

        pErrb = new teestream::teebuf(pLogBuf, orig_cerrb);
        std::cerr.rdbuf(pErrb);
    }

    // info lines
    {
        std::string status(aef_git_status);
        bool bdirty = status.contains('M') || status.contains('d');
        std::string dirty = bdirty ? "dirty" : "clean";
        std::cout << "Stark-diagonalizer version compiled on " << __DATE__ << " "
            << __TIME__ << ", git commit " << aef_git_commit << std::endl;
        std::cout << "Git status is " << dirty << " string {" << status << "}" << std::endl;
        std::cout << fmt::format("Start time is {}", start_time) << std::endl;
        std::cout << fmt::format("Eigen will use {} threads", Eigen::nbThreads()) << std::endl;
    }
#ifdef _OPENMP
    std::cout << "Reconfiguring openmp to use the correct number of threads (the number of physical cores)." << std::endl;
    int num_physical_cores = get_num_cores();
    omp_set_num_threads(num_physical_cores);
    Eigen::setNbThreads(num_physical_cores);
    std::cout << fmt::format("OpenMP/Eigen will use {} threads", num_physical_cores) << std::endl;
#endif

#ifndef DONT_USE_CUDA
    constexpr bool diag_use_cuda = true;
    std::cout << "Initializing matrix backend" << std::endl;
    aef::ResultCode rc = aef::matrix::init(aef::matrix::BackendType::NvidiaCuda, argc, argv);
    if (!aef::succeeded(rc)) {
        std::cout << fmt::format("Initializing matrix backend failed with error {} = 0x{:x}", static_cast<int32_t>(rc), static_cast<uint32_t>(rc));
    }
    std::cout << "Successfully initialized CUDA" << std::endl;
#else
    constexpr bool diag_use_cuda = false;
    aef::matrix::init(aef::matrix::BackendType::EigenCPU, argc, argv);
#endif

    std::cout << "Constructing HyperfineCalculator with nmax = " << param_nmax << std::endl;
    std::cout << fmt::format("nmax is {}, E_z is {} MHz/D, K is {} MHz ({})",
        param_nmax, calc_E_z, 0, "disabled") << std::endl;
    HyperfineCalculator calc(param_nmax, calc_E_z);
    std::cout << "Finished making HyperfineCalculator with nmax = " << param_nmax << std::endl;
#ifndef DONT_USE_CUDA
    Eigen::MatrixXcd vals;
    vals.resize(calc.nBasisElts, calc.nBasisElts);
    vals.setZero();

    std::cout << fmt::format(
        "Setting up matrix backend device-side buffers with nRows={} after creating molecular system",
        calc.nBasisElts) << std::endl;
    aef::matrix::set_max_size(calc.nBasisElts);
#endif


    prev_time = log_time_at_point("Starting matrix element calculations", start_time, prev_time);
    calc.calculate_matrix_elts();
    //calc.diagonalize_H(diag_use_cuda);

    prev_time = log_time_at_point("Finished matrix elt calcs", start_time, prev_time);
    
    // set H_tot = H_stk
    calc.H_tot = calc.H_stk;
    calc.diagonalize_H(diag_use_cuda);
    prev_time = log_time_at_point("Diagonalized Stark Potential", start_time, prev_time);

    if (param_nmax >= 20) {
        calc.save_matrix_elts(dpath / "matrix.dat");
        prev_time = log_time_at_point("Saved hyprfinecalculator", start_time, prev_time);
    }
#if 0
    // diagonalize H_stk
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
    solver.compute(calc.H_stk);
    calc.Es = solver.eigenvalues();
    calc.Vs = solver.eigenvectors();
#endif

    // output state coeffs
    {
        fs::path csvpath = odir / fmt::format("{}.csv", calc_E_z);
        std::ofstream out(csvpath);
        // write out header line
        out << "Eidx";
        for (int kdx = 0; kdx < calc.nBasisElts; kdx++) {
            out << fmt::format(",Re(<j_{}|E_n>),Im(<j_{}|E_n>", kdx, kdx);
        }
        out << std::endl;

        for (int Edx = 0; Edx < calc.nBasisElts; Edx++) {
            dump_state(calc, Edx, out);
        }
        out.close();
    }
    // output MDA info
    {
        std::ofstream oState(odir / fmt::format("state_ifo_{}.csv", calc_E_z));
        output_state_info(oState, calc
#ifndef DONT_USE_CUDA
            , vals
#endif
        );
        oState.close();
    }

    // Output low_state_dumper style files
    {
        fs::path sdir = odir / "state_coeffs";
        fs::create_directories(sdir);
        fs::path csvpath = sdir / fmt::format("{}.csv", calc_E_z_V_per_cm);
        std::ofstream out(csvpath);
        // write out header line
        out << "Eidx";
        for (int kdx = 0; kdx < calc.nBasisElts; kdx++) {
            out << fmt::format(",Re(<j_{}|E_n>),Im(<j_{}|E_n>", kdx, kdx);
        }
        out << std::endl;
        // dump +Z-oriented f=0/f=1 singlet-triplet
        for (int idx = bidx_posz; idx < bidx_posz + num_singlet_triplet; idx++) {
            dump_state(calc, idx, out);
        }
        // dump -Z-oriented f=0/f=1 singlet-triplet
        for (int idx = bidx_negz; idx < bidx_negz + num_singlet_triplet; idx++) {
            dump_state(calc, idx, out);
        }

        // dump "+-XY oriented" f=0/f=1 singlet-triplets (note: these mix --> not actually oriented)
        for (int idx = bidx_pmxy; idx < bidx_pmxy + 4 * num_singlet_triplet; idx++) {
            dump_state(calc, idx, out);
        }
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

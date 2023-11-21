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

namespace fs = std::filesystem;
using namespace std::chrono;

#include "../AeF-hyperfine-structure.inl"

time_point<system_clock> log_time_at_point(
    const char* desc,
    time_point<system_clock>& start,
    time_point<system_clock>& prev) {
    using namespace std::chrono;
    time_point<system_clock> curr_time = system_clock::now();
    std::string stime = fmt::format("{0:%F}-{0:%H%M}{0:%S}", curr_time);
    using d_seconds = duration<double>;
    // get seconds elapsed since start_time
    auto st_diff = curr_time - start;
    double st_sec_count = d_seconds(st_diff).count();
    // calculate seconds elapsed since prev_time
    auto pv_diff = curr_time - prev;
    double pv_sec_count = d_seconds(pv_diff).count();

    std::string lstr = fmt::format(
        "{}: have taken {} seconds since last, {} seconds since start (current time is {})", desc,
        pv_sec_count, st_sec_count, stime);
    std::cout << lstr << std::endl;
    *(std::addressof(prev)) = curr_time;
    return curr_time;
}

/// <summary>
/// Dumps a row from a given matrix.  This
/// </summary>
/// <param name="Vs">Matrix to</param>
/// <param name="state"></param>
/// <param name="out"></param>
void dump_row(Eigen::MatrixXcd& Vs, int state, std::ostream& out) {
    out << state;
    for (int kdx = 0; kdx < Vs.cols(); kdx++) {
        out << "," << std::real(Vs(state, kdx)) << "," << std::imag(Vs(state, kdx));
    }
    out << std::endl;
}

/// <summary>
/// Performs the quantum q*(q+1) squaring for angular momenta
/// </summary>
/// <param name="q"> q </param>
/// <returns>q*(q+1)</returns>
template <class T> T qsq(T q) {
    return q * (q + (T)1);
}


/// <summary>
/// Inverts the quantum q*(q+1) squaring for angular momenta
/// </summary>
/// <param name="expect_qsq">the expectation value of &lt;q(q+1)&gt; </param>
/// <returns>the effective expectation of q</returns>
template <class T> T invert_qsq(T expect_qsq) {
    return (std::sqrt(4 * expect_qsq + 1.0) - 1.0) / 2.0;
}

/// <summary>
/// Calculates the expectation values of the basis operators.  The squared
/// </summary>
/// <param name="calc">HyperfineCalculator: contains operator matrix elements
/// and states</param> <param name="E_idx">the index of Energy level to
/// calculate with</param> <returns></returns>
j_basis_vec expectation_values_jsq(HyperfineCalculator& calc, int32_t E_idx) {
    //
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    j_basis_vec out;
#ifdef _WIN32
    SecureZeroMemory((void*)&out, sizeof(j_basis_vec));
#else
    //explicit_bzero((void*)&out, sizeof(j_basis_vec));
    memset((void*)&out, 0, sizeof(j_basis_vec));
#endif
    double prob_tot = 0;
    for (int32_t kidx = 0; kidx < calc.nBasisElts; kidx++) {
        const double prob = std::norm(state_vec[kidx]);

        if (prob < std::numeric_limits<double>::epsilon()) {
            continue;
        }

        prob_tot += prob;
        j_basis_vec bs_ket = calc.basis[kidx];
        // note that angular momenta
        out.n += prob * bs_ket.n * (bs_ket.n + 1);
        out.j += prob * bs_ket.j * (bs_ket.j + 1);
        out.f += prob * bs_ket.f * (bs_ket.f + 1);
        out.m_f += prob * bs_ket.m_f;
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }

    out.n = invert_qsq(out.n / prob_tot);
    out.j = invert_qsq(out.j / prob_tot);
    out.f = invert_qsq(out.f / prob_tot);
    out.m_f /= prob_tot;

    return out;
}

double expect_parity(HyperfineCalculator& calc, int32_t E_idx) {
    double ex_parity = 0.0;
    double prob_tot = 0.0;
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    for (int32_t kidx = 0; kidx < calc.nBasisElts; kidx++) {
        const double prob = std::norm(state_vec[kidx]);

        if (prob < std::numeric_limits<double>::epsilon()) {
            continue;
        }

        prob_tot += prob;
        j_basis_vec bs_ket = calc.basis[kidx];
        ex_parity += prob * std::pow(-1, bs_ket.n);
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }
    return ex_parity / prob_tot;
}

/// <summary>
/// Outputs some key information about each energy eigenstate including:
/// * the MDA expectation value (re then IM)
/// * the expectation values of n, j, f, and m_f
/// </summary>
/// <param name="output">the stream to output to</param>
/// <param name="calc"></param>
void output_state_info(std::ostream& output, HyperfineCalculator& calc) {
    output << "Index n, Energy (MHz), Re(<n|dx|n>), Re(<n|dy|n>), Re(<n|dz|n>), "
        "Im(<n|dx|n>), Im(<n|dy|n>), Im(<n|dz|n>), "
        "Re(<n|n|n>), Re(<n|j|n>), Re(<n|f|n>), Re(<n|m_f|n>),"
        "Im(<n|n|n>), Im(<n|j|n>), Im(<n|f|n>), Im(<n|m_f|n>),"
        "<n|(-1)^n|n>"
        << std::endl;

    for (size_t n = 0; n < calc.nBasisElts; n++) {
        auto e_n = calc.Vs.col(n);
        // molecular dipole vector in spherical tensor form
        dcomplex d10 = expectation_value(e_n, calc.d10);
        dcomplex d11 = expectation_value(e_n, calc.d11);
        dcomplex d1t = expectation_value(e_n, calc.d1t);

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
    /// 50 kV/cm = 25170 MHz/D is the field strength used to calculate H_stk.
    /// DO NOT CHANGE THIS.
    /// This is a standardized value used across the AeF-hyperfine-structure codebase
    /// </summary>
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
    tee::teebuf* pBuf, * pOutb, * pErrb;

    std::streambuf* pLogBuf;
    std::streambuf* orig_coutb = std::cout.rdbuf();
    std::streambuf* orig_cerrb = std::cerr.rdbuf();

    {
        auto fpath = dpath / "out.log";
        oLog = std::ofstream(fpath, std::ios::trunc | std::ios::out);
        // test if logging to debug window enabled
        if (enable_debug_log) {
            // if enabled, tee to both the log file and the debug window
            pBuf = new tee::teebuf(oLog.rdbuf(), pDebug->rdbuf());
            pLogBuf = pBuf;
        } else {
            // otherwise just output to the log file
            pBuf = nullptr;
            pLogBuf = oLog.rdbuf();
        }
        pOutb = new tee::teebuf(pLogBuf, orig_coutb);
        std::cout.rdbuf(pOutb);

        pErrb = new tee::teebuf(pLogBuf, orig_cerrb);
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




    HyperfineCalculator calc(param_nmax, calc_E_z);

    prev_time = log_time_at_point("Starting matrix element calculations", start_time, prev_time);
    calc.calculate_matrix_elts();
    calc.diagonalize_H();
    if (param_nmax >= 20)
        calc.save_matrix_elts(dpath / "matrix.dat");

    prev_time = log_time_at_point("Finished matrix elt calcs", start_time, prev_time);

    // set H_tot = H_stk
    calc.H_tot = calc.H_stk;
    calc.diagonalize_H();
#if 0
    // diagonalize H_stk
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
    solver.compute(calc.H_stk);
    calc.Es = solver.eigenvalues();
    calc.Vs = solver.eigenvectors();
#endif

    fs::path csvpath = odir / fmt::format("{}.csv", calc_E_z);
    std::ofstream out(csvpath);
    // write out header line
    out << "Eidx";
    for (int kdx = 0; kdx < calc.nBasisElts; kdx++) {
        out << fmt::format(",Re(<j_{}|E_n>),Im(<j_{}|E_n>", kdx, kdx);
    }
    out << std::endl;

    for (int Edx = 0; Edx < calc.nBasisElts; Edx++) {
        dump_row(calc.Vs, Edx, out);
    }
    out.close();
    std::ofstream oState(odir / fmt::format("state_ifo_{}.csv", calc_E_z));
    output_state_info(oState, calc);
    oState.close();
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

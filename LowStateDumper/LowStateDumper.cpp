// LowStateDumper.cpp : This file contains the 'main' function. Program execution begins and ends there.
/*
    LowStateDumper/LowStateDumper.cpp -- this program dumps some of the 
    low-lying states of a loaded HyperfineCalculator.  For the case of 138BaF
    when the devonshire potential is enabled, these states correspond to the
    lowest lying positive-z and negative-z singlet-triplet group. 

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
#include <fmt.hpp>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <cxxopts.hpp>

using namespace std::chrono;
namespace fs = std::filesystem;

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

void dump_state(HyperfineCalculator& calc, int state, std::ostream &out) {
    out << state;
    // BUGFIX: states are column vectors not row vectors
    Eigen::VectorXcd state_vec = calc.Vs.col(state);
    for (int kdx = 0; kdx < calc.nBasisElts; kdx++) {
        dcomplex ampl = state_vec[kdx];
        out << "," << std::real(ampl) << "," << std::imag(ampl);
    }
    out << std::endl;
}

#include "../AeF-hyperfine-structure.inl"

int main(int argc, char **argv) {
    //auto dpath = fs::path("output");
    std::chrono::time_point<std::chrono::system_clock> start_time =
        std::chrono::system_clock::now();
    std::string stime = fmt::format("{0:%F}-{0:%H%M}{0:%S}", start_time);
    std::chrono::time_point<std::chrono::system_clock> prev_time = start_time;
    //dpath /= stime;
    //fs::create_directories(dpath);

    bool enable_debug_log = false;
    std::string loadname = "matrix.dat";
    bool print_extras = true;
    bool output_Es = true;
    std::vector<double> E_zs;
    bool k_specified = false;
    double K = 0;
    bool k_enabled = false;

    /// <summary>
    /// The value of E_z used to calculate H_stk is 50 kV/cm = 25170 MHz/D.
    /// This isn't saved to file but all correctly-calculated aefdats use this
    /// value.
    /// </summary>
    constexpr double calc_E_z = unit_conversion::MHz_D_per_V_cm * 50 * 1000;
    
    cxxopts::Options options("low-state-dumper", "Program to dump the low-lying +-Z oriented states of interest");

    options.add_options()
        ("h,help", "Print usage")
        ("e,Ez", "Electric field values to consider in V/cm (can use multiple times)", cxxopts::value<std::vector<double>>()) // allow multiple Ez
        ("d,enable_debug", "Enable debug mode", cxxopts::value<bool>()->default_value("false"))
        ("print_extras", "Print extra information", cxxopts::value<bool>()->default_value("true"))
        ("l,load", "Load molecular system operators from file", cxxopts::value<std::string>())
        ("o,output", "Set output directory.  Default is coeffs", cxxopts::value<std::string>())
        ("k", "K value to use in Kelvin -- version 2 AeFDat files don't have this in them, although it is contained in the log]", cxxopts::value<double>());
    
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("enable_debug")) {
        enable_debug_log = result["enable_debug"].as<bool>();
    }

    if (result.count("load")) {
        loadname = result["load"].as<std::string>();
    }

    if (result.count("Ez")) {
        using std::begin, std::end;
        auto vals = result["Ez"].as<std::vector<double>>();
        E_zs.insert(E_zs.end(), begin(vals), end(vals));
    } else {
        E_zs.push_back(500);
    }

    if (result.count("k")) {
        k_specified = true;
        K = result["k"].as<double>() * unit_conversion::MHz_per_Kelvin;
        k_enabled = K != 0;
    }
    
    fs::path p = fs::absolute(loadname);
    fs::path dir = p.parent_path();
    
    fs::path odir = dir / "state_coeffs";
    fs::create_directories(odir);

    // set up logging
    std::ofstream oLog(odir / "low_state_dumper.log", std::ios::trunc | std::ios::out);
    aef::LogRedirector lredir(oLog, false, true);

    {
        std::string status(aef_git_status);
        bool bdirty = status.contains('M') || status.contains('d');
        std::string dirty = bdirty ? "dirty" : "clean";
        std::cout << "AeF Hyperfine Structure LowStateDumper, version compiled on " << __DATE__ << " "
            << __TIME__ << ", git commit " << aef_git_commit << std::endl;
        std::cout << "Git status is " << dirty << " string {" << status << "}" << std::endl;
        std::cout << fmt::format("Start time is {}", start_time) << std::endl;
        std::cout << fmt::format("Eigen will use {} threads", Eigen::nbThreads()) << std::endl;
    }

    prev_time = log_time_at_point("[Low state dumper] Initializing", start_time, prev_time);

#ifndef DONT_USE_CUDA
    constexpr bool diag_use_cuda = true;
    std::cout << "Initializing CUDA" << std::endl;
    aef::init_cuda(argc, (const char**)argv);
    std::cout << "Successfully initialized CUDA" << std::endl;
#else
    constexpr bool diag_use_cuda = false;
#endif

#ifdef _OPENMP
    std::cout << "Reconfiguring openmp to use the correct number of threads (the number of physical cores)." << std::endl;
    int num_physical_cores = get_num_cores();
    omp_set_num_threads(num_physical_cores);
    Eigen::setNbThreads(num_physical_cores);
    std::cout << fmt::format("OpenMP/Eigen will use {} threads", num_physical_cores) << std::endl;
#endif

    std::cout << fmt::format("[Low State dumper] attempting to load matrix elements from {}", p.string()) << std::endl;
    
    HyperfineCalculator calc(4, 0, K);
    bool success = calc.load_matrix_elts(loadname);

    if (!success) {
        // load failed, attempt to diagnose why & exit
        prev_time = log_time_at_point("[Low state dumper] loading matrix elements failed, see above for details", start_time, prev_time);
        bool exists = fs::exists(loadname);
        bool isfile = fs::is_regular_file(loadname);
        
        if (!exists) {
            std::cout << std::format("Hint: \"{}\" appears not to exist", loadname) << std::endl;
        } else if (!isfile) { // not relevant if file doesn't exist
            std::cout << std::format("Hint: \"{}\" appears to exist but not be a regular file", loadname) << std::endl;
        } else if (!have_read_perms(loadname)) {
            std::cout << std::format("Hint: \"{}\" appears to be unreadable but exists and is a regular file.", loadname) << std::endl;
        }  else {
            std::cout << std::format("Hint: \"{}\" appears to exist, be a regular file, and be readable, so it's probably corrupt", loadname) << std::endl;
        }
        std::exit(255);
    }

#ifndef DONT_USE_CUDA
        std::cout << fmt::format(
            "Setting up CUDA device-side buffers with nRows={} after loading matrix elements",
            calc.nBasisElts) << std::endl;
        aef::cuda_resize(calc.nBasisElts);
        std::cout << "Finished CUDA device-side buffer setup" << std::endl;
#endif


    // note: we actually don't need the actual value of K here, since H_dev is stored pre-multiplied by K
    // we only need to know whether K is zero (in-vacuum) or not (in-matrix).
    if (k_specified) {
        // Starting with aefdat_version::rawmat_okq_params, aefdat files save E_z and K, thus specifying 
        if (calc.load_version >= aefdat_version::rawmat_okq_params) {
            std::cout << fmt::format("Attempt") << std::endl;
            exit(253);
        }
    }
    std::cout << fmt::format("K_enabled = {}", calc.enableDev ? "true" : "false") << std::endl;
    // XXX this is only needed because AeFDat used to not save K
    if (!k_specified && calc.load_version <= aefdat_version::rawmat_okq_params) {
        // attempt to locate out.log
        fs::path olog_path = dir / "out.log";
        
        if (!fs::is_regular_file(olog_path)) {
            auto str = std::format("K neither specified nor present in AeFDat matrix file, and was unable locate log path {}",
                olog_path.string());
            std::cout << str << std::endl;
            exit(254);
        }
        // TODO parse log file???

    }

    std::cout << "does H_tot commute with F_z? " << aef::commutes(calc.H_tot, calc.F_z) << std::endl;
    std::cout << "does H_stk commute with F_z? " << aef::commutes(calc.H_stk, calc.F_z) << std::endl;
    std::cout << "does H_hfs commute with F_z? " << aef::commutes(calc.H_hfs, calc.F_z) << std::endl;

    // load succeded
    prev_time = log_time_at_point("[Low state dumper] Finished loading matrix elements.", start_time, prev_time);

    // give info --> note: enableDev and K weren't saved in older formats 
    std::cout << fmt::format("Molecular System information: nmax = {}", calc.nmax) << std::endl;//, calc.enableDev);
    std::cout << fmt::format("[] Starting to dump states, {} field values", E_zs.size()) << std::endl;

    for (auto Ez : E_zs) {
        const double Ez_mhz = Ez * unit_conversion::MHz_D_per_V_cm;
        std::string desc = fmt::format("[Low state dumper] starting calculations for electric field {}", Ez);
        prev_time = log_time_at_point(desc.c_str(), start_time, prev_time);
        std::cout << "Diagonalizing matrix" << std::endl;
        const double scale = Ez_mhz / calc_E_z;
        std::cout << fmt::format("H_tot rows={}, cols={}", calc.H_tot.rows(), calc.H_tot.cols()) << std::endl;
        std::cout << fmt::format("H_rot rows={}, cols={}", calc.H_rot.toDenseMatrix().rows(), calc.H_rot.toDenseMatrix().cols()) << std::endl;
        std::cout << fmt::format("H_hfs rows={}, cols={}", calc.H_hfs.rows(), calc.H_hfs.cols()) << std::endl;
        std::cout << fmt::format("H_stk rows={}, cols={}", calc.H_stk.rows(), calc.H_stk.cols()) << std::endl;
        calc.H_tot = calc.H_rot.toDenseMatrix() + calc.H_hfs + scale * calc.H_stk + calc.H_dev;

        calc.diagonalize_H(diag_use_cuda);

        fs::path csvpath = odir / fmt::format("{}.csv", Ez);
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
        desc = fmt::format("[Low state dumper] Finished calculations for electric field {}", Ez);
        prev_time = log_time_at_point(desc.c_str(), start_time, prev_time);
        out.close();
    }
    lredir.touch();
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

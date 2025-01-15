// PerturbationAnalyzer.cpp : This file contains the 'main' function. Program execution begins and ends there.
// AeF-hyperfine-structure.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
// This code implements the "main program" of the aef-hyperfine-structure toolkit.
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
#include <pch.h>
#include <aef/aef.h>
#include <aef/debug_stream.h>
#include <aef/matrix_utils.h>
#include <aef/teestream.hpp>
#include <aef/aef_run.h>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fmt.hpp>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <cxxopts.hpp>
#include <aef/quantum.h>
#include <aef/aef_run.h>
#include <aef/operators/operators.h>


using namespace std::chrono;
namespace fs = std::filesystem;
namespace hfs_constants = baf_constants;
using aef::log_time_at_point;
using namespace aef::quantum;

#include "../AeF-hyperfine-structure.inl"



int main(int argc, char **argv) {

    constexpr double calc_E_z = unit_conversion::MHz_D_per_V_cm * 50 * 1000;

    std::chrono::time_point<std::chrono::system_clock> start_time =
        std::chrono::system_clock::now();
    std::string stime = fmt::format("{0:%F}-{0:%H%M}{0:%S}", start_time);
    std::chrono::time_point<std::chrono::system_clock> prev_time = start_time;
    //fs::create_directories(dpath);

    int param_nmax = 20;
    bool enable_debug_log = false;
    bool load_from_file = false;
    std::string loadname = "";
    bool print_extras = true;
    bool output_Es = true;
    size_t nStarkIterations = 101;
    double min_E_z = 0;
    double max_E_z = calc_E_z;
    fs::path dpath("output");

    // todo parse args
    // args should include: E_max, nmax, enable_debug_log
    cxxopts::Options options("aef-hyperfine-structure", "Program to simulate the hyperfine structure of"
        "diatomic Alkaline - monoflouride molecules");
    options.add_options()
        ("h,help", "Print usage")
        ("e,E_max", "Maximum electric field [V/cm]", cxxopts::value<double>())
        ("E_min", "Minimum electric field [V/cm]", cxxopts::value<double>())
        ("n,n_max", "Maximum n level to include", cxxopts::value<int>())
        ("d,enable_debug", "Enable debug mode", cxxopts::value<bool>()->default_value("false"))
        ("print_extras", "Print extra information", cxxopts::value<bool>()->default_value("true"))
        ("l,load", "Load molecular system operators from file", cxxopts::value<std::string>())
        ("t,stark_iterations", "Number of iterations to perform the stark loop for", cxxopts::value<size_t>());

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

    if (result.count("stark_iterations")) {
        nStarkIterations = result["stark_iterations"].as<size_t>();
    }

    if (result.count("E_max")) {
        max_E_z = result["E_max"].as<double>();
    }

    if (result.count("E_min")) {
        min_E_z = result["E_min"].as<double>();
    }

    if (!load_from_file) {
        std::clog << "[PerturbationAnalyzer] Error: must load from file" << std::endl;
    }
    fs::path runpath = aef::get_aef_run_path(fs::absolute(loadname));
    dpath = runpath / "ptfw";
    fs::create_directories(dpath);

    // create info log
    std::ofstream oLog(dpath / "perturbation_analyzer_out.log", std::ios::trunc | std::ios::out);
    aef::LogRedirector lredir(oLog, enable_debug_log, true);
    // info lines
    {
        std::string status(aef_git_status);
        bool bdirty = status.contains('M') || status.contains('d');
        std::string dirty = bdirty ? "dirty" : "clean";
        std::cout << "AeF Hyperfine Structure perturbation analyzer, version compiled on " << __DATE__ << " "
            << __TIME__ << ", git commit " << aef_git_commit << std::endl;
        std::cout << "Git status is " << dirty << " string {" << status << "}" << std::endl;
        std::cout << fmt::format("Start time is {}", start_time) << std::endl;
        std::cout << fmt::format("Eigen will use {} threads", Eigen::nbThreads()) << std::endl;
    }

    // log arguments
    {
        std::cout << "Arguments: [";
        for (int i = 0; i < argc; i++) {
            std::cout << fmt::format(" {{{}}}", argv[i]);
        }
        std::cout << "]" << std::endl;
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



    init_rng();

    aef::aef_run run(runpath);

    HyperfineCalculator calc;
    bool succ = calc.load_matrix_elts(run.get_matrix_path());

    if (!succ) {
        // TODO error
        std::string str_pth = run.get_matrix_path().generic_string();
        std::cerr << fmt::format("[Perturbation Analyzer] Loading matrix file {} from run {} failed!",
            str_pth, run.get_run_name()) << std::endl;
        std::abort();
        aef::unreachable();
    }

#ifndef DONT_USE_CUDA
    Eigen::MatrixXcd vals;
    vals.resize(calc.nBasisElts, calc.nBasisElts);
    vals.setZero();

    std::cout << fmt::format(
        "Setting up matrix backend device-side buffers with nRows={} after creating molecular system",
        calc.nBasisElts) << std::endl;
    aef::matrix::set_max_size(calc.nBasisElts);
#endif

    rc = aef::ResultCode::Success;

    // make bigmatrix
    prev_time = log_time_at_point("Creating perturbation theory framework", start_time, prev_time);
    aef::operators::PerturbationFramework pfw(&calc);
    prev_time = log_time_at_point("Constructing eEDM-like operator", start_time, prev_time);
    pfw.addOperator("eEDM", new aef::operators::eEDMOperator());
    prev_time = log_time_at_point("Constructing 19F NSM-like operator", start_time, prev_time);
    pfw.addOperator("NSM", new aef::operators::NSMOperator());
    prev_time = log_time_at_point("Constructing Z-axis Stark operator", start_time, prev_time);
    pfw.addOperator("StarkZ", new aef::operators::StarkOperator({0.0,0.0,1.0}));
    prev_time = log_time_at_point("Constructing Z-xis Zeeman operator", start_time, prev_time);
    pfw.addOperator("ZeemanZ", new aef::operators::ZeemanOperator({0, 0, 1.0}));

    prev_time = log_time_at_point("Evaluating operators", start_time, prev_time);
    pfw.evaluate();
    prev_time = log_time_at_point("Finished evaluating operators", start_time, prev_time);

    /// Start delta
    prev_time = log_time_at_point("Starting perturbative energy shift calculations", start_time, prev_time);
    // delta E vector from eEDM-like operator
    Eigen::VectorXcd dEs_eEDM;
    prev_time = log_time_at_point("1st ord PT eEDM start", start_time, prev_time);
    rc = pfw.delta_E_lo("eEDM", dEs_eEDM);
    if (!aef::succeeded(rc)) {
        // error
        std::clog << fmt::format("delta-E eEDM calc failed {}", (int)rc);
    }
    prev_time = log_time_at_point("1st ord PT eEDM done", start_time, prev_time);

    // delta E vector from Fluorine-19 nuclear schiff moment-like operator
    prev_time = log_time_at_point("1st ord PT NSM start", start_time, prev_time);
    Eigen::VectorXcd dEs_f_nsm; 
    rc = pfw.delta_E_lo("NSM", dEs_f_nsm);
    if (!aef::succeeded(rc)) {
        // error
        std::clog << fmt::format("delta-E light NSM calc failed {}", (int)rc);
    }
    prev_time = log_time_at_point("1st ord PT NSM done", start_time, prev_time);

    // delta-E vector from Z-axis Zeeman shift
    prev_time = log_time_at_point("1st ord PT ZeemanZ start", start_time, prev_time);
    Eigen::VectorXcd dEs_zeez;
    rc = pfw.delta_E_lo("ZeemanZ", dEs_zeez);
    if (!aef::succeeded(rc)) {
        // error
        std::clog << fmt::format("delta-E Z-axis Zeeman calc failed {}", (int)rc);
    }
    prev_time = log_time_at_point("1st ord PT ZeemanZ done", start_time, prev_time);
    prev_time = log_time_at_point("Perturbative energy shift calculations complete", start_time, prev_time);

    // file output
    std::cout << "Energy Eigenstate Index\tDelta E eEDM (MHz)\tDelta E 19F NSM (MHz)\tDelta E Zeeman Z (MHz)" << std::endl;
    for (int idx = 0; idx < calc.nBasisElts; idx++) {
        std::cout << fmt::format("{}\t{}\t{}", idx, dEs_eEDM(idx), dEs_f_nsm(idx), dEs_zeez(idx)) << std::endl;
    }

    std::ofstream out(dpath / "cpv_energies.tsv");
    out << "Energy Eigenstate Index\tDelta E eEDM (MHz)\tDelta E 19F NSM (MHz)\tDelta E Z-axis Zeeman (MHz)"
        "\tImaginary Part of dE_EDM(MHz)\tImaginary Part of dE_19F_NSM(MHz)\tImagninary Part of dE_ZeeZ (MHz)" << std::endl;
    for (int idx = 0; idx < calc.nBasisElts; idx++) {
        dcomplex dE_EDM = dEs_eEDM(idx);
        dcomplex dE_f_NSM = dEs_f_nsm(idx);
        dcomplex dE_zeez = dEs_zeez(idx);
        out << fmt::format("{}\t{}\t{}\t{}\t{}", idx, std::real(dE_EDM), std::real(dE_f_NSM), std::real(dE_zeez),
            std::imag(dE_EDM), std::imag(dE_f_NSM), std::imag(dE_zeez)) << std::endl;
    }
    out.close();
    prev_time = log_time_at_point("File writes complete4", start_time, prev_time);
//    std::cout << fmt::format("Energy vector {}", es) << std::endl;
    aef::matrix::shutdown();
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

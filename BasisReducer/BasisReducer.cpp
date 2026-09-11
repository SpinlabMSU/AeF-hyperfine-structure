// BasisReducer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
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
#include <system_error>
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
#include <aef/MolecularSystem.h>


using namespace std::chrono;
namespace fs = std::filesystem;
using aef::log_time_at_point;
using namespace aef::quantum;

#include "../AeF-hyperfine-structure.inl"


using namespace std::chrono;
namespace fs = std::filesystem;
using system_time = std::chrono::time_point<std::chrono::system_clock>;

//int32_t closest_approx(Eigen)

namespace {
    enum class orient_choice {
        INVALID = 0,
        FORCE_DISABLE = 1,
        AUTO = 2,
        FORCE_ENABLE = 3
    };
    std::istream& operator>>(std::istream& in, orient_choice& c) {
        int a;
        in >> a;
        c = (orient_choice)(a + 1);
        return in;
    }
};
/*
enum qoperator:: {

    reduceAndOutputOperator(sys, sys.H_tot, (char*)"h_tot", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperatorRe(sys, sys.H_tot, (char*)"h_tot", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.H_stk, (char*)"h_stk", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.H_dev, (char*)"h_dev", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.H_hfs, (char*)"h_hfs", work, rbasis_size, dpath, start_time, &prev_time);
    vals = sys.H_rot.toDenseMatrix();
    reduceAndOutputOperator(sys, vals, (char*)"h_rot", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.d10, (char*)"O_d10", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.d11, (char*)"O_d11", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.d1t, (char*)"O_d1t", work, rbasis_size, dpath, start_time, &prev_time);
};
*/
struct reduced_system {
    int num;
    Eigen::MatrixXcd* h_tot;
};

void reduceMatrix(Eigen::MatrixXcd& opReducedOut, Eigen::MatrixXcd& opJfBasis, aef::MolecularSystem &sys, int size, Eigen::MatrixXcd *work) {
    //Eigen::MatrixXcd reducedVs = sys.Vs(Eigen::seq())
    // eigen is column major
    // for now just reduce at the end -- todo 
    opReducedOut.resize(size, size); opReducedOut.setZero();
    aef::matrix::group_action(*work, sys.Vs, opJfBasis);
    opReducedOut = (*work)(Eigen::seq(0, size-1), Eigen::seq(0, size-1));
}

void reduceAndOutputOperator(aef::MolecularSystem& sys, Eigen::MatrixXcd& op, char* fnam, Eigen::MatrixXcd& work, int size, fs::path p, 
    system_time start_time, system_time *prev_time) {
    Eigen::MatrixXcd opReduced;

    *prev_time = log_time_at_point(fmt::format("Reducing operator {} to size {}", fnam, size).c_str(), start_time, *prev_time);
    reduceMatrix(opReduced, op, sys, size, &work);

    *prev_time = log_time_at_point("Reducing finished, now writing reduced operator", start_time, *prev_time);
    std::string spath = fmt::format("{}.csv", fnam);
    std::ofstream out(p / spath);
    out << fnam;
    for (int jdx = 0; jdx < size; jdx++) {
        // write out column index
        out << fmt::format(", {}", jdx);
    }
    out << std::endl;
    // 
    for (int idx = 0; idx < size; idx++) {
        // write out row index
        out << fmt::format("{}", idx);
        for (int jdx = 0; jdx < size; jdx++) {
            // write out element -- remember that Eigen is column major
            auto elt = opReduced(jdx, idx);
            out << fmt::format(", ({}+i*{})", std::real(elt), std::imag(elt));
        }
        out << std::endl;
    }
    *prev_time = log_time_at_point(fmt::format("Done with operator {}", fnam).c_str(), start_time, *prev_time);
}

void reduceAndOutputOperatorRe(aef::MolecularSystem& sys, Eigen::MatrixXcd& op, char* fnam, Eigen::MatrixXcd& work, int size, fs::path p,
    system_time start_time, system_time* prev_time) {
    Eigen::MatrixXcd opReduced;

    *prev_time = log_time_at_point(fmt::format("Reducing operator {} to size {}", fnam, size).c_str(), start_time, *prev_time);
    reduceMatrix(opReduced, op, sys, size, &work);

    *prev_time = log_time_at_point("Reducing finished, now writing reduced operator real part", start_time, *prev_time);
    std::string spath = fmt::format("{}_re.csv", fnam);
    std::ofstream out(p / spath);
    out << fnam;
    for (int jdx = 0; jdx < size; jdx++) {
        // write out column index
        out << fmt::format(", {}", jdx);
    }
    out << std::endl;
    // 
    for (int idx = 0; idx < size; idx++) {
        // write out row index
        out << fmt::format("{}", idx);
        for (int jdx = 0; jdx < size; jdx++) {
            // write out element -- remember that Eigen is column major
            auto elt = opReduced(jdx, idx);
            out << fmt::format(", {}", std::real(elt));
        }
        out << std::endl;
    }
    *prev_time = log_time_at_point(fmt::format("Done with operator {}", fnam).c_str(), start_time, *prev_time);
}

using aef::quantum::transition_type;

struct tsn_op_file {
    std::string path_add;
    std::ofstream stream;
};



void reduceAndOutputTransitionOperator(aef::MolecularSystem& sys,
    Eigen::MatrixXcd& op, transition_type typ, int order, 
    char* fnam, Eigen::MatrixXcd& work, int size, fs::path dir,
    system_time start_time, system_time* prev_time) {
    Eigen::MatrixXcd opReduced;
    *prev_time = log_time_at_point(fmt::format("Reducing operator {} to size {} with magsq", fnam, size).c_str(), start_time, *prev_time);
    reduceMatrix(opReduced, op, sys, size, &work);

    *prev_time = log_time_at_point("Reducing finished, now writing reduced operator and magsq", start_time, *prev_time);
    

    std::string spath = fmt::format("{}.csv", fnam);
    std::ofstream out(dir / spath);

    std::string spath2 = fmt::format("{}_magsq.csv", fnam); // op mag sq
    std::ofstream out2(dir / spath2);
    
    std::string spathA = fmt::format("{}_A.csv", fnam); // einstein A coeff
    std::ofstream outA(dir / spathA);

    std::string spathMag = fmt::format("{}_mag.csv", fnam); // operator magnitude
    std::ofstream outMag(dir / spathMag);

    std::string spathf = fmt::format("{}_f.csv", fnam); // einstein A coeff
    std::ofstream outf(dir / spathf);

    std::vector<std::ofstream*> files = {
        &out, &out2, &outA, &outMag, &outf
    };

    for (int idx = 0; idx < files.size(); idx++) {
        std::ofstream& of = *files[idx];
        of << fmt::format("{}", fnam);
        for (int jdx = 0; jdx < size; jdx++) {
            of << fmt::format(", {}", jdx);
        }
        of << std::endl;
    }
    // 
    for (int idx = 0; idx < size; idx++) {
        // write out row index
        for (int jdx = 0; jdx < files.size(); jdx++) {
            std::ofstream& of = *files[jdx];
            of << fmt::format("{}", idx);
        }
        for (int jdx = 0; jdx < size; jdx++) {
            // write out element -- remember that Eigen is column major
            auto elt = opReduced(jdx, idx);
            out << fmt::format(", ({}+i*{})", std::real(elt), std::imag(elt));
            out2 << fmt::format(", {}", std::norm(elt));
            outMag << fmt::format(", {}", std::abs(elt));

            const double S = std::norm(elt);
            const double E_i = std::real(sys.Es(idx));
            const double E_j = std::real(sys.Es(jdx));
            const double dE = std::abs(std::real(E_i - E_j));
            assert(!isnan(dE));
            aef::quantum::transition_information tsn(typ, order, dE,elt);
            double A = tsn.calc_A(); //aef::quantum::calculate_transition_rate(typ, order, dE, elt);
            double f = tsn.calc_f();
            const char* tsn_type = (typ == transition_type::E) ? "E" : "M";
            constexpr double thresh = 1E-4;
            if (std::abs(elt) > thresh){
                std::string strI = "";
                std::string strJ = "";

                if (!sys.enableDev) {
                    strI = fmt::format("{}", expectation_values_jsq(sys, idx));
                    strJ = fmt::format("{}", expectation_values_jsq(sys, jdx));
                }

                std::cout << fmt::format("{8} - {0}{1} transition {2} {3}\n"
                    "\t{2}{10}: E{2} = {4} MHz\n"
                    "\t{3}{11}: E{3} = {5} MHz\n"
                    "\tdE = {6} MHz,\tS={9} D^2\n"
                    "\tA = {7} Hz, f = {12}\n", tsn_type, order, idx, jdx, E_i, E_j, dE, A, fnam, S, strI, strJ, f);
            }
            outA << fmt::format(", {}", A);
            outf << fmt::format(", {}", f);
        }
        for (int jdx = 0; jdx < files.size(); jdx++) {
            std::ofstream& of = *files[jdx];
            of << std::endl;
        }
    }
    for (int jdx = 0; jdx < files.size(); jdx++) {
        std::ofstream& of = *files[jdx];
        of.close();
    }
    *prev_time = log_time_at_point(fmt::format("Done with operator {}", fnam).c_str(), start_time, *prev_time);
}


int main(int argc, char **argv) {
    constexpr std::string_view progname("BasisReducer");
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
    double E_z = calc_E_z;
    double E_z_V_cm = calc_E_z / unit_conversion::MHz_D_per_V_cm;
    bool E_z_specified = false;
    fs::path dpath("output");
    double B_z_Gauss = 0;
    double B_x_Gauss = 0;
    double B_y_Gauss = 0;
    bool B_specified = false;
    int rbasis_size = 48;

    orient_choice o_choice = orient_choice::AUTO;

    // todo parse args
    // args should include: E_max, nmax, enable_debug_log
    cxxopts::Options options("basis-reducer", "Program to calculate perturbative corrections to"
        " the hyperfine structure of diatomic Alkaline - monofluoride molecules");
    options.add_options()
        ("h,help", "Print usage")
        ("e,Ez", "Electric field for PT calculations [V/cm]", cxxopts::value<double>())
        ("E_min", "Minimum electric field [V/cm]", cxxopts::value<double>())
        ("b,Bz", "Magnetic field along the Z-axis for PT calculations [G]", cxxopts::value<double>())
        ("Bx", "Magnetic field along the X-axis for PT calculations [G]", cxxopts::value<double>())
        ("By", "Magnetic field along the Y-axis for PT calculations [G]", cxxopts::value<double>())
        ("n,n_max", "Maximum n level to include", cxxopts::value<int>())
        ("r,reduced_size", "Size of the reduced-basis system", cxxopts::value<int>())
        ("d,enable_debug", "Enable debug mode", cxxopts::value<bool>()->default_value("false"))
        ("print_extras", "Print extra information", cxxopts::value<bool>()->default_value("true"))
        ("l,load", "Load molecular system operators from file", cxxopts::value<std::string>())
        ("t,stark_iterations", "Number of iterations to perform the stark loop for", cxxopts::value<size_t>())
        ("O,orientation_diagonalizer", "Force Enable (2)/Force Disable(0)/Auto(1) Orientation Diagonalizer", cxxopts::value<orient_choice>());

    options.allow_unrecognised_options();
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

    if (result.count("Ez")) {
        E_z_V_cm = result["Ez"].as<double>();
        E_z_specified = true;
    }

    if (result.count("E_min")) {
        min_E_z = result["E_min"].as<double>();
    }

    if (result.count("Bz")) {
        B_z_Gauss = result["Bz"].as<double>();
        B_specified = true;
    }

    if (result.count("Bx")) {
        B_x_Gauss = result["Bx"].as<double>();
        B_specified = true;
    }

    if (result.count("By")) {
        B_y_Gauss = result["By"].as<double>();
        B_specified = true;
    }

    if (result.count("orientation_diagonalizer")) {
        o_choice = result["orientation_diagonalizer"].as<orient_choice>();
    }

    if (result.count("reduced_size")) {
        rbasis_size = result["reduced_size"].as<int>();
    }

    if (!load_from_file) {
        std::clog << "[" << progname <<"] Error: must load from file" << std::endl;
        exit(1);
    }

    if (rbasis_size <= 0) {
        std::clog << 
            fmt::format("[{}] Error: reduced basis size {} is less than or equal to zero",
                progname, rbasis_size) 
            << std::endl;
        exit(1);
    }

    std::error_code ec;
    aef::aef_run run(fs::absolute(loadname));
    fs::path runpath = run.get_run_path();
    dpath = (runpath / "reduced") / fmt::format("{}", rbasis_size);// / fmt::format("{}", );

    if (B_specified) {
        E_z = E_z_V_cm * unit_conversion::MHz_D_per_V_cm;
        dpath /= fmt::format("E_{}_B_{}_{}_{}", E_z_V_cm, B_x_Gauss, B_y_Gauss, B_z_Gauss);
    } else if (E_z_specified) {
        E_z = E_z_V_cm * unit_conversion::MHz_D_per_V_cm;
        dpath /= fmt::format("{}", E_z_V_cm);
    }

    std::cout << fmt::format("[{}] Using output directory {}", progname, dpath.generic_string()) << std::endl;
    if (!fs::exists(dpath)) {
        fs::create_directories(dpath, ec);
        if (ec) {
            std::clog << fmt::format("[{}] Unable to create output directory, error category: {}, code: {}, message: {}",
                progname, ec.category().name(), ec.value(), ec.message()) << std::endl;
            exit(2);
        }
    } else if (!fs::is_directory(dpath)) {
        //
        std::clog << fmt::format("Error: output path {} exists but is not a directory!", dpath.string()) << std::endl;
    }

    // create info log
    std::ofstream oLog(dpath / "basis_reducer.log", std::ios::trunc | std::ios::out);
    aef::LogRedirector lredir(oLog, enable_debug_log, true);
    // info lines
    {
        std::string status(aef_git_status);
        bool bdirty = status.contains('M') || status.contains('d');
        std::string dirty = bdirty ? "dirty" : "clean";
        std::cout << "AeF Hyperfine Structure basis reducer, version compiled on " << __DATE__ << " "
            << __TIME__ << ", git commit " << aef_git_commit << std::endl;
        std::cout << "Git status is " << dirty << " string {" << status << "}" << std::endl;
        std::cout << fmt::format("Start time is {}", start_time) << std::endl;
        std::cout << fmt::format("Eigen will use {} threads", Eigen::nbThreads()) << std::endl;
        std::string Ez_spec = E_z_specified ? "" : " not";
        std::cout << fmt::format("E_z has{} been specified, E_z = {} MHz/D = {} V/cm", Ez_spec, E_z, E_z / unit_conversion::MHz_D_per_V_cm) << std::endl;
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
    prev_time = log_time_at_point("Initializing matrix Backend", start_time, prev_time);
#ifndef DONT_USE_CUDA
    constexpr bool diag_use_cuda = true;
    std::cout << fmt::format("{} Initializing matrix backend", progname) << std::endl;
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
    std::cout << "Successfully initialized RNG" << std::endl;

    //aef::aef_run run(runpath);
    prev_time = log_time_at_point("Loading molecular system", start_time, prev_time);
    aef::MolecularSystem sys;
    {
        auto mpath = run.get_run_path() / "molsys.dat";
        std::string str_pth = mpath.generic_string();
        std::cout << fmt::format("[{}] Loading matrix file from {}", progname, str_pth) << std::endl;
        rc = sys.load(mpath); //run.get_matrix_path());

        if (aef::failed(rc)) {
            // TODO error
            std::string str_pth = run.get_matrix_path().generic_string();
            std::cerr << fmt::format("[{}] Loading matrix file {} from run {} failed!",
                progname, str_pth, run.get_run_name()) << std::endl;
            std::abort();
            aef::unreachable();
        }
    }


    Eigen::MatrixXcd vals, work;
    vals.resize(sys.nBasisElts, sys.nBasisElts);
    vals.setZero();
    work.resize(sys.nBasisElts, sys.nBasisElts);
    work.setZero();

    prev_time = log_time_at_point("Matrix backend setup", start_time, prev_time);
    std::cout << fmt::format(
        "Setting up matrix backend device-side buffers with nRows={} after creating molecular system",
        sys.nBasisElts) << std::endl;
    aef::matrix::set_max_size(sys.nBasisElts);
    prev_time = log_time_at_point("Recalculating H_tot with specified E_z", start_time, prev_time);

    // 
    Eigen::MatrixXcd H_zeeman;
    H_zeeman.resize(sys.nBasisElts, sys.nBasisElts);
    H_zeeman.setZero();
    {
        aef::operators::ZeemanOperator zop(sys, { 0, 0, 1 });
        zop.fillMatrix(H_zeeman);
    }

    Eigen::MatrixXcd Dev_orient_Diagonalizer;
    Dev_orient_Diagonalizer.setZero();
    if (sys.enableDev) {
        Dev_orient_Diagonalizer = aef::orient_diag::makeOrientationDiagonalizer(sys);
    }

    if (B_specified || E_z_specified) {
        prev_time = log_time_at_point("Recalculating H_tot with specified E_z and B_z", start_time, prev_time);
        const double scale = E_z / calc_E_z;
        sys.H_tot = sys.H_rot.toDenseMatrix() + sys.H_hfs + scale * sys.H_stk + sys.H_dev + H_zeeman;
        prev_time = log_time_at_point("Finished recalculating H_tot, now diagonalizing", start_time, prev_time);
        if (sys.enableDev) {
            (void)aef::orient_diag::diagonalize(sys, Dev_orient_Diagonalizer, &vals);
        } else {
            sys.diagonalize();
        }
        prev_time = log_time_at_point("Diagonalization complete", start_time, prev_time);
    }

    {
        // output unit file
        std::ofstream out(dpath / "units.txt", std::ios::trunc | std::ios::out);
        out << "UNIT LIST\n-----------------" << std::endl;
        out << "Hamiltonians and Energies: MHz (i.e actually the _linear_ frequency)" << std::endl;
        out << "Molecular Dipolar Axis d1p: dimensionless" << std::endl;
        out << "E1 transition dipoles E1_d1p: Debye" << std::endl;
        out << "M1 transition dipoles M1_d1p: (MHz/T)" << std::endl;
        out << "E1 transition dipole Magnitude: Debye" << std::endl;
        out << "M1 transition dipole Magnitude: (MHz/T)" << std::endl;
        out << "E1 transition dipole Squared Magnitude: Debye ^ 2" << std::endl;
        out << "M1 transition dipole Squared Magnitude: (MHz/T) ^ 2" << std::endl;
        out << "E1 transition Einstein A: Hz" << std::endl;
        out << "M1 transition Einstein A: Hz" << std::endl;
        out << "E1 transition Einstein B: Hz" << std::endl;
        out << "M1 transition Einstein B: Hz" << std::endl;
        out << "E1 transition f/Oscillator Strength: Dimensionless" << std::endl;
        out << "M1 transition f/Oscillator Strength: Dimensionless" << std::endl;
        out.close();
    }

    rc = aef::ResultCode::Success;

    
    reduceAndOutputOperator(sys, sys.H_tot, (char*)"h_tot", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperatorRe(sys, sys.H_tot, (char*)"h_tot", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.H_stk, (char*)"h_stk", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.H_dev, (char*)"h_dev", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys, sys.H_hfs, (char*)"h_hfs", work, rbasis_size, dpath, start_time, &prev_time);
    vals = sys.H_rot.toDenseMatrix();
    reduceAndOutputOperator(sys,      vals, (char*)"h_rot", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys,   sys.d10, (char*)"O_d10", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys,   sys.d11, (char*)"O_d11", work, rbasis_size, dpath, start_time, &prev_time);
    reduceAndOutputOperator(sys,   sys.d1t, (char*)"O_d1t", work, rbasis_size, dpath, start_time, &prev_time);
    
    Eigen::MatrixXcd vals2, vals3;
    vals2.resize(sys.nBasisElts, sys.nBasisElts); vals2.setZero();
    vals3.resize(sys.nBasisElts, sys.nBasisElts); vals3.setZero();
    {
        // for E1 transitions
        sys.get_calc()->calculate_mol_EDM(vals, vals2, vals3);
        using transition_type::E;
        reduceAndOutputTransitionOperator(sys, vals , E, 1, (char*)"E1_d10", work, rbasis_size, dpath, start_time, &prev_time);
        reduceAndOutputTransitionOperator(sys, vals2, E, 1, (char*)"E1_d11", work, rbasis_size, dpath, start_time, &prev_time);
        reduceAndOutputTransitionOperator(sys, vals3, E, 1, (char*)"E1_d1t", work, rbasis_size, dpath, start_time, &prev_time);
    }

    {
        // for M1 transitions
        sys.get_calc()->calculate_mol_MDM(vals, vals2, vals3);
        using transition_type::M;
        static_assert(transition_type::E != transition_type::M);
        reduceAndOutputTransitionOperator(sys, vals , M, 1, (char*)"M1_d10", work, rbasis_size, dpath, start_time, &prev_time);
        reduceAndOutputTransitionOperator(sys, vals2, M, 1, (char*)"M1_d11", work, rbasis_size, dpath, start_time, &prev_time);
        reduceAndOutputTransitionOperator(sys, vals3, M, 1, (char*)"M1_d1t", work, rbasis_size, dpath, start_time, &prev_time);
    }
    return 0;
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

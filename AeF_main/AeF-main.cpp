// AeF-main.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
// This code implements the new "main program" of the aef-hyperfine-structure toolkit.
// The modification is that this code now uses the aef::MolecularSystem code
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
#include "AeF-main.h"

#include "AeF-main.inl"
// int32_t closest_approx(Eigen)



////// Stateful dot-product state tracking: spinlab memebers see elog:EDM3/680
// note: these will be empty when tracking is disabled
bool _do_tracking;
Eigen::MatrixXcd prev_Vs_H; // store previous eigenvectors in hermition conjugate form
Eigen::VectorXcd prev_Es;
std::vector<size_t> sdx_from_edx_map; // sdx_from_edx_map[energy eigenstate idx] == state idx
std::vector<size_t> edx_from_sdx_map; // edx_from_sdx_map[state idx] = energy eigenstate idx

void init_state_tracking(aef::MolecularSystem& sys, double stk_scale=1.0) {
    _do_tracking = true;
    // construct as
    sys.H_tot = sys.H_rot.toDenseMatrix() + sys.H_hfs + stk_scale * sys.H_stk + sys.H_dev;
    sys.diagonalize();
    prev_Vs_H = sys.Vs.adjoint();
    prev_Es = sys.Es;
    const size_t nBasisElts = sys.nBasisElts;
    // initial mapping will be the identity mapping
    sdx_from_edx_map.resize(sys.nBasisElts);
    edx_from_sdx_map.resize(sys.nBasisElts);
    for (size_t idx = 0; idx < nBasisElts; idx++) {
        sdx_from_edx_map[idx] = idx;
        edx_from_sdx_map[idx] = idx;
    }
}

void disable_tracking() {
    _do_tracking = false;
    prev_Es.resize(0);
    prev_Vs_H.resize(0, 0);
    sdx_from_edx_map.clear();
    edx_from_sdx_map.clear();
}

void update_tracking(aef::MolecularSystem& sys, Eigen::MatrixXcd &scratch) {
    // tracks new against old
    const size_t nBasisElts = sys.nBasisElts;
    Eigen::MatrixXcd& Vs = sys.Vs;
    Eigen::VectorXcd& Es = sys.Es;

    // Probability matrix is an matrix W_{ij} = || < E_old_i | E_new_j > ||^2
    aef::matrix::multiply(prev_Vs_H, Vs, scratch);
    Eigen::MatrixXcd& probs = scratch;
    probs = probs.array().abs2();

    std::unordered_set<size_t> used_values;
    double worst_max_prob = std::numeric_limits<double>::infinity();
    Eigen::Index worst_edx = -1, worst_sdx = -1;

    /// This implementation is completely unoptimized
    for (size_t sdx = 0; sdx < nBasisElts; sdx++) {
        // evaluate "best" match
        double max_prob = -std::numeric_limits<double>::infinity();
        Eigen::Index max_edx = -1;
        for (size_t edx = 0; edx < nBasisElts; edx++) {
            double prob = std::real(probs(sdx, edx)); // discard imaginary part, it should be zero
            if (prob > max_prob && !used_values.contains(edx)) {
                max_edx = edx;
                max_prob = prob;
            }
        }
        assert((max_edx != -1) && "Unable to find best match (this should be mathematically impossible)");
        edx_from_sdx_map[sdx] = max_edx;
        sdx_from_edx_map[max_edx] = sdx;
        used_values.insert(max_edx);

        if (worst_max_prob >= max_prob) {
            worst_max_prob = max_prob;
            worst_edx = max_edx;
            worst_sdx = sdx;
        }
    }

    std::cout << fmt::format("Worst max probability was {} at edx = {}, sdx = {}", worst_max_prob, worst_edx,
                worst_sdx) << std::endl;

    // finish by setting the previous values to the current values so we can update the current values
    prev_Vs_H = sys.Vs.adjoint();
    prev_Es = sys.Es;
}

/// <summary>
/// This writes a CSV file that holds both 
/// </summary>
void output_tracking_info(std::ostream &out) {
    out << "Index,State Idx From Energy E-state Index,Energy E-state Index from State Index" << std::endl;
    for (size_t idx = 0; idx < edx_from_sdx_map.size(); idx++) {
        out << fmt::format("{},{},{}", idx, sdx_from_edx_map[idx], edx_from_sdx_map[idx]) << std::endl;
    }
    out.flush();
}

/// <summary>
/// Calculates the expectation values of an energy eigenstate
/// </summary>
/// <param name="calc">HyperfineCalculator: contains operator matrix elements
/// and states</param> <param name="E_idx">the index of Energy level to
/// calculate with</param> <returns></returns>
aef::universal_diatomic_basis_vec expectation_values(aef::MolecularSystem& calc, int32_t E_idx) {
    //
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    aef::universal_diatomic_basis_vec out;
    aef::IMolecularCalculator *mcalc = calc.get_calc();
#ifdef _WIN32
    SecureZeroMemory((void*)&out, sizeof(aef::universal_diatomic_basis_vec));
#else
    // explicit_bzero isn't neccessary and reduces portability
    //explicit_bzero((void*)&out, sizeof(aef::universal_diatomic_basis_vec));
    memset((void*)&out, 0, sizeof(aef::universal_diatomic_basis_vec));
#endif
    out = mcalc->get_basis_ket(0);
    out.r = out.n = out.j = out.f_1 = out.f = out.m_f = 0;
    double prob_tot = 0;
#ifdef MATRIX_ELT_DEBUG
    double expect_mf = -0.1;
    bool set_mf = false;
#endif
    for (int32_t kidx = 0; kidx < calc.nBasisElts; kidx++) {
        const double prob = std::norm(state_vec[kidx]);

        if (prob < std::numeric_limits<double>::epsilon()) {
            continue;
        }

        prob_tot += prob;
        aef::universal_diatomic_basis_vec bs_ket = mcalc->get_basis_ket(kidx);

#ifdef MATRIX_ELT_DEBUG
        if (prob > 1) {
            DebugBreak();
        }

        if (set_mf && expect_mf != calc.basis[kidx].m_f && prob > 0) {
            DebugBreak();
            throw(false);
        }

        if (prob > 0 && !set_mf) {
            expect_mf = bs_ket.m_f;
            set_mf = true;
        }
#endif
        out.r += prob * bs_ket.r;
        out.n += prob * bs_ket.n;
        out.j += prob * bs_ket.j;
        out.f_1 += prob * bs_ket.f_1;
        out.f += prob * bs_ket.f;
        out.m_f += prob * bs_ket.m_f;
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }

    out.r /= prob_tot;
    out.n /= prob_tot;
    out.j /= prob_tot;
    out.f_1 /= prob_tot;
    out.f /= prob_tot;
    out.m_f /= prob_tot;

    return out;
}

static inline double diff_states(aef::universal_diatomic_basis_vec v1, aef::universal_diatomic_basis_vec v2) {
    // parameters need to be tweaked for best results
    constexpr double cn = 3.5;
    constexpr double cj = 0.25; // 1.0;
    constexpr double cf1 = 0.5;
    constexpr double cf = 2.5;
    constexpr double cm = 5.0; // 100.0;

    double dn = (v1.n - v2.n);
    double dj = (v1.j - v2.j);
    double df1 = (v1.f_1 - v2.f_1);
    double df = (v1.f - v2.f);
    double dm = (v1.m_f - v2.m_f);
    return cn * dn * dn + cj * dj * dj + cf1 * df1 * df1 + cf * df * df + cm * dm * dm;
}

/// <summary>
/// Finds the diagonlized eignestate "closest" to the basis state with the provided index
/// "ket_idx"
/// </summary>
/// <param name="calc">System basis and operator elements</param>
/// <param name="ket_idx">the index of the desired basis state</param>
/// <param name="exclude_Eidx">(optional) an energy eigenstate to "exclude" from
/// being the closest state.  Intended to fix the </param> <returns></returns>
int32_t closest_state(aef::MolecularSystem& calc, int32_t ket_idx,
    int32_t exclude_Eidx = -1) {
    int32_t closest_idx = -1;
//#define USE_EXPECTATION_VALUES
#define USE_TRIVIAL
#define USE_STATE_TRACKING
#if defined(USE_STATE_TRACKING)
    if (_do_tracking) {
        closest_idx = sdx_from_edx_map[ket_idx];
    } else {
        closest_idx = ket_idx; // tracking disable -- just do trivial
    }
#elif defined(USE_EXPECTATION_VALUES
    // strategy 1: look for the energy eigenstate whose expectation values most closely match the target state
    double chisq = (double)std::numeric_limits<double>::infinity();
    aef::universal_diatomic_basis_vec ket = calc.get_calc()->get_basis_ket(ket_idx);

    for (int32_t Eidx = 0; Eidx < calc.nBasisElts; Eidx++) {
        aef::universal_diatomic_basis_vec expect_ket = expectation_values(calc, Eidx);
        double localx2 = diff_states(ket, expect_ket);

        if (localx2 < chisq && Eidx != exclude_Eidx) {
            chisq = localx2;
            closest_idx = Eidx;
        }
    }
#elif defined(USE_TRIVIAL)
    // strategy 2: just assume the indices stay the same
    closest_idx = ket_idx;
#else
    // strategy 3: assume the best state will have the largest probability of ket_idx
    Eigen::VectorXcd ket_coeffs = calc.Vs.row(ket_idx); // vector of <E_idx|ket_idx>
    double max_magsq = -1;
    for (int32_t Eidx = 0; Eidx < calc.nBasisElts; Eidx++) {
        double magsq = std::norm(ket_coeffs[Eidx]);
        if (magsq > max_magsq) {
            max_magsq = magsq;
            closest_idx = Eidx;
        }
    }
#endif
    return closest_idx;
}

/// <summary>
/// Outputs some key information about each energy eigenstate including:
/// * the MDA expectation value (re then IM)
/// * the expectation values of n, j, f, and m_f
/// </summary>
/// <param name="output">the stream to output to</param>
/// <param name="calc"></param>
void output_state_info(std::ostream& output, aef::MolecularSystem& calc
#ifndef DONT_USE_CUDA
    , Eigen::MatrixXcd &vals
#endif
) {
    output << "Index n, Energy (MHz), Re(<n|dx|n>), Re(<n|dy|n>), Re(<n|dz|n>), "
        "Im(<n|dx|n>), Im(<n|dy|n>), Im(<n|dz|n>), "
        "Re(<n|n|n>), Re(<n|j|n>), Re(<n|f_1|n>), Re(<n|f|n>), Re(<n|m_f|n>),"
        "Im(<n|n|n>), Im(<n|j|n>), Im(<n|f_1|n>), Im(<n|f|n>), Im(<n|m_f|n>),"
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
        aef::universal_diatomic_basis_vec v = expectation_values_jsq(calc, n);

        // output
        auto mda_ifo =
            fmt::format("{}, {}, {}, {}, {}, {}, {}, {}", n, std::real(calc.Es[n]),
                std::real(dx), std::real(dy), std::real(dz), std::real(dx),
                std::imag(dy), std::imag(dz));
        auto re_njf1fmf = fmt::format("{},{},{}, {},{}", std::real(v.n), std::real(v.j), std::real(v.f_1), std::real(v.f), std::real(v.m_f));
        auto im_njf1fmf = fmt::format("{},{},{}, {},{}", std::imag(v.n), std::imag(v.j), std::imag(v.f_1), std::imag(v.f), std::imag(v.m_f));
        output << mda_ifo << ", " << re_njf1fmf << ", " << im_njf1fmf << "," << expect_parity(calc, n) << std::endl;
    }
    output.flush();
}

int main(int argc, char **argv) {
    ///////////////////// constants
    /// <summary>
    /// 50 kV/cm = 25170 MHz/D is the field strength used to calculate H_stk.
    /// DO NOT CHANGE THIS.
    /// Instead use the E_max argument
    /// </summary>
    constexpr double calc_E_z = unit_conversion::MHz_D_per_V_cm * 50 * 1000;

#ifndef _MAKEFILE_PROVIDES_DEVFLAG
#define USE_DEVONSHIRE
#endif

#ifdef USE_DEVONSHIRE
    /// <summary>
    /// Devonshire coupling constant: 100 Kelvin --> 2.083 THz
    /// </summary>
    constexpr double K = 100.0 * unit_conversion::MHz_per_Kelvin;
    constexpr const char* devstatus = "enabled";
#else
    /// <summary>
    /// Disable Devonshire potential: 0 MHz --> disabled
    /// </summary>
    constexpr double K = 0;
    constexpr const char* devstatus = "disabled";
#endif

    ///////////////////////// main code

    // output file --> automatically make output based on current datetime
    auto dpath = fs::path("output");
    std::chrono::time_point<std::chrono::system_clock> start_time =
        std::chrono::system_clock::now();
    std::string stime = fmt::format("{0:%F}-{0:%H%M}{0:%S}", start_time);
    std::chrono::time_point<std::chrono::system_clock> prev_time = start_time;
    dpath /= stime;
    // fs::create_directories(dpath); // don't perform the directory creation until after parsing arguments
    // this prevents running --help from creating empty useless directories

    int param_nmax = 20;
    bool enable_debug_log = false;
    bool load_from_file = false;
    std::string loadname = "";
    bool print_extras = true;
    bool output_Es = true;
    bool force_save = false;
    size_t nStarkIterations = 101;
    double min_E_z = 0;
    double max_E_z = calc_E_z / unit_conversion::MHz_D_per_V_cm; // units of max_E_z are V/cm
    std::string mol_calc_type = aef::RaFMolecularCalculator::calc_type_str;
    bool do_tracking = true;
    bool round_files = true;

    // todo parse args
    // args should include: E_max, nmax, enable_debug_log
    cxxopts::Options options("aef-hyperfine-structure", "Program to calculate the hyperfine structure of"
        " diatomic Alkaline-monofluoride molecules");
    options.add_options()
        ("h,help", "Print usage")
        ("e,E_max", "Maximum electric field [V/cm]", cxxopts::value<double>())
        ("Z,E_min", "Minimum electric field [V/cm]", cxxopts::value<double>())
        ("n,n_max", "Maximum n level to include", cxxopts::value<int>())
        ("d,enable_debug", "Enable debug mode", cxxopts::value<bool>()->default_value("false"))
        ("print_extras", "Print extra information", cxxopts::value<bool>()->default_value("true"))
        ("l,load", "Load molecular system operators from file", cxxopts::value<std::string>())
        ("t,stark_iterations", "Number of iterations to perform the stark loop for", cxxopts::value<size_t>())
        ("s,sys", "Molecular system type to use", cxxopts::value<std::string>())
        ("S,force-save", "Force the Molecular System to always be saved", cxxopts::value<bool>()->default_value("false"))
        ("do_tracking", "Do state tracking", cxxopts::value<bool>()->default_value("true"))
        ("round_files", "Round the electric field to the nearest integer in the info_Ez_{}.csv files", cxxopts::value<bool>()->default_value("true"));

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

    if (result.count("sys")) {
        mol_calc_type = result["sys"].as<std::string>();
    }

    if (result.count("force-save")) {
        force_save = result["force-save"].as<bool>();
    }

    if (result.count("do_tracking")) {
        do_tracking = result["do_tracking"].as<bool>();
    }

    if (result.count("round_files")) {
        round_files = result["round_files"].as<bool>();
    }

    // Create output directory and info log now that arguments have been parsed
    fs::create_directories(dpath);
    std::ofstream oLog(dpath / "out.log", std::ios::trunc | std::ios::out);
    aef::LogRedirector lredir(oLog, enable_debug_log, true);
    // info lines
    {
        std::string status(aef_git_status);
        bool bdirty = status.contains('M') || status.contains('d');
        std::string dirty = bdirty ? "dirty" : "clean";
        std::cout << "AeF-Hyperfine-Structure main spectrum calculation program (MolecularSystem enhanced), version compiled on " << __DATE__ << " "
            << __TIME__ << ", git commit " << aef_git_commit << ", main program file " __FILE__ << std::endl;
        std::cout << "Git status is " << dirty << " string {" << status << "}" << std::endl;
        std::cout << fmt::format("Start time is {}", start_time) << std::endl;
        std::cout << fmt::format("Eigen will use {} threads", Eigen::nbThreads()) << std::endl;
    }

    // log arguments
    {
        std::cout << "Arguments: [";
        for (int i = 0; i < argc; i++){
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

    init_rng();
#ifndef DONT_USE_CUDA
    constexpr bool diag_use_cuda = true;
    std::cout << "Initializing matrix backend" << std::endl;
    aef::matrix::init(aef::matrix::BackendType::NvidiaCuda, argc, argv);
    std::cout << "Successfully initialized CUDA" << std::endl;
#else
    constexpr bool diag_use_cuda = false;
    aef::matrix::init(aef::matrix::BackendType::EigenCPU, argc, argv);
#endif

    // maximum value of the n quantum number.  There will be 8*(nmax**2) states
    int nmax = param_nmax;
    aef::IMolecularCalculator* pCalc = nullptr;

    constexpr bool enable_dynamic_mol_calc = true;
    if(enable_dynamic_mol_calc) {
        pCalc = aef::IMolecularCalculator::makeCalculatorOfType(mol_calc_type); //new aef::RaFMolecularCalculator(nmax);
        if (!pCalc) {
            std::cerr << fmt::format("Error: calculator type \"{}\" does not exist", mol_calc_type) << std::endl;
            exit(-98);
        } else {
            std::clog << fmt::format("Using calculator type \"{}\", actual type name \"{}\"", mol_calc_type, pCalc->get_calc_type());
        }
        pCalc->set_nmax(nmax);
    } else {
        pCalc = new aef::BaFMolecularCalculator(nmax);
    }

    aef::MolecularSystem sys(pCalc, nmax, calc_E_z, K);
    {
        std::string track_status = do_tracking ? "enabled" : "disabled";
        std::cout << fmt::format("nmax is {}, E_z is {} MHz/D, K is {} MHz ({}), calculator is {}, tracking is {}",
            nmax, calc_E_z, K, devstatus, pCalc->get_calc_type(), track_status) << std::endl;
    }
    if (load_from_file) {
        std::string logstr = fmt::format("Loading matrix elements from {}", loadname);
        prev_time = log_time_at_point(logstr.c_str(), start_time, prev_time);
        aef::ResultCode res = sys.load(loadname);

        if (aef::failed(res)) {
            std::cout << "couldn't load " << loadname << " error code " << static_cast<uint32_t>(res) << std::endl;
            exit(-1);
        }
        logstr = fmt::format("Finished loading matrix elements from {}", loadname);
        prev_time = log_time_at_point(logstr.c_str(), start_time, prev_time);
#ifndef DONT_USE_CUDA
        std::cout << fmt::format(
            "Setting up CUDA device-side buffers with nRows={} after loading matrix elements",
            sys.nBasisElts) << std::endl;
        aef::matrix::set_max_size(sys.nBasisElts);
        std::cout << "Finished CUDA device-side buffer setup" << std::endl;
#endif
    } else {
        // not loading from file --> calculate
#ifndef DONT_USE_CUDA
        std::cout << fmt::format(
            "Setting up backend device-side buffers with nRows={} before matrix element calculations",
            sys.nBasisElts) << std::endl;
        aef::matrix::set_max_size(sys.nBasisElts);
        std::cout << "Finished backend device-side buffer setup" << std::endl;
#endif
        {
            prev_time = log_time_at_point("Starting matrix element calculations", start_time, prev_time);
            sys.calculate_matrix_elts();
            auto prev_2 = log_time_at_point("[Not updating global previous time] finished actual matrix element calculations", start_time, prev_time);
                sys.diagonalize();
                prev_2 = log_time_at_point("[Not updating global previous time] finished actual matrix element calculations", start_time, prev_2);
                if (force_save || nmax >= 20) {
                    sys.save(dpath / "molsys.dat");
                    prev_2 = log_time_at_point("[Not updating previous time globally] finished saving molecular system", start_time, prev_2);
                }
            prev_time = log_time_at_point("Finished matrix elt calcs", start_time, prev_time);
        }
    }

    if (print_extras) {
        Eigen::VectorXcd Es = sys.Es;
        std::cout << "----------- Stark-Shifted -----------" << std::endl;
        std::cout << "Level, Energy (MHz)" << std::endl;
        double EPrev = 0;
        for (int i = 0; i < sys.nBasisElts; i++) {
            double dE = std::real(Es[i]) - EPrev;
            std::cout << i << ", " << std::real(Es[i]) << "MHz, DeltaE = " << dE
                << " MHz, " << expectation_values(sys, i).ket_string()
                << std::endl;
            EPrev = std::real(Es[i]);
        }
        std::cout << std::endl << std::endl;

        std::cout << "----------- NO STARK -----------" << std::endl;
        sys.H_tot = sys.H_rot.toDenseMatrix() + sys.H_hfs; // -= sys.H_stk;
        sys.diagonalize();

        Es = sys.Es;
        EPrev = 0;
        std::cout << "Level, ket, Energy (MHz)" << std::endl;
        for (int i = 0; i < sys.nBasisElts; i++) {
            double dE = std::real(Es[i]) - EPrev;
            std::cout << i << ", " << std::real(Es[i]) << "MHz, DeltaE = " << dE
                << " MHz, " << expectation_values(sys, i).ket_string() << std::endl;
            EPrev = std::real(Es[i]);
            // std::cout << i << ", " << calc.basis[i] << ", " << Es[i] << std::endl;
        }
    } else {
        sys.H_tot -= sys.H_stk;
        sys.diagonalize();
    }
#ifndef DONT_USE_CUDA
    Eigen::MatrixXcd vals;
    vals.resize(sys.nBasisElts, sys.nBasisElts);
    vals.setZero();
#endif

    // create output file
    auto fpath = dpath / "stark_shift_gnd.csv";
    std::ofstream oStk(fpath, std::ios::trunc | std::ios::out);

    auto epath = dpath / "stark_spectrum.csv";
    std::ofstream oEs (epath, std::ios::trunc | std::ios::out);
    oEs << "idx,E-field (V/cm)";
    for (size_t idx = 0; idx < sys.nBasisElts; idx++) {
        oEs << fmt::format(",E{}", idx);
    }
    oEs << std::endl;

    aef::universal_diatomic_basis_vec gnd = pCalc->get_basis_ket(0);


    std::vector<aef::universal_diatomic_basis_vec> lowest_states = pCalc->get_lowest_states();
    const int nLowestStates = lowest_states.size();
    std::vector<int> lowest_idxs(nLowestStates);

    for (int sdx = 0; sdx < nLowestStates; sdx++) {
        lowest_idxs[sdx] = pCalc->get_index(lowest_states[sdx]);
    }

    // oStk << "E-field (V/cm), Stark-shifted Energy of " << gnd.ket_string() << "(MHz)";
    assert(sys.H_tot.rows() == sys.H_stk.rows());

    // output stark.csv header line + lowest states
    oStk << "E-field (V/cm), dE_gnd";// << ", dE_f1t, dE_f10, dE_f11" << std::endl;
    for (int sdx = 1; sdx < nLowestStates; sdx++) {
        std::cout << fmt::format("Lowest state #{}: ket #{}, {}", sdx, lowest_idxs[sdx], lowest_states[sdx]) << std::endl;
        oStk << ", dE_" << sdx;
    }
    oStk << std::endl;


#define USE_REAL_Es
#ifdef USE_REAL_Es
    typedef double etype;
#define EVAL(val) std::real(val) 
#else
    typedef dcomplex etype;
#define EVAL(val) val
#endif

    Eigen::MatrixXcd lowest_energies(nStarkIterations, nLowestStates);

#ifndef USE_DEVONSHIRE
    // note: devonshire potential doesn't conserve m_f
    for (int idx = 0; idx < sys.nBasisElts; idx++) {
        aef::universal_diatomic_basis_vec v1 = pCalc->get_basis_ket(idx);//calc.basis[idx];
        for (int jdx = 0; jdx < sys.nBasisElts; jdx++) {
            aef::universal_diatomic_basis_vec v2 = pCalc->get_basis_ket(jdx);
            double prob = std::norm(sys.H_tot(idx, jdx));
            if (v1.m_f != v2.m_f && prob > 0) {
                DebugBreak();
                std::string ostr =
                    fmt::format("ERROR v1 = {}, v2 = {}, prob = {}", v1, v2, prob);
                std::cout << ostr << std::endl;
                std::cerr << ostr << std::endl;
                assert(!(v1.m_f != v2.m_f && prob > 0));
            }
        }
    }
#else
#endif // !USE_DEVONSHIRE
    // directory to put devonshire info
    auto devpath = dpath / "devonshire_info";
    fs::create_directories(devpath);
    std::cout << "does H_tot commute with d10? " << aef::matrix::commutes(sys.H_tot, sys.d10, &vals) << std::endl;
    std::cout << "does H_tot commute with d11? " << aef::matrix::commutes(sys.H_tot, sys.d11, &vals) << std::endl;
    std::cout << "does H_tot commute with d1t? " << aef::matrix::commutes(sys.H_tot, sys.d1t, &vals) << std::endl;
    std::cout << std::endl;

    std::cout << "Is d10  all zero " << sys.d10.isZero(1E-6) << std::endl;
    std::cout << "Is d11  all zero " << sys.d11.isZero(1E-6) << std::endl;
    std::cout << "Is d1t  all zero " << sys.d1t.isZero(1E-6) << std::endl;
    std::cout << "Is Hdev all zero " << sys.H_dev.isZero(1E-6) << std::endl;

    // Initialize state tracking
    auto track_dir_path = dpath / "tracking_info";
    if (do_tracking) {
        prev_time = log_time_at_point("Initializing State Tracking", start_time, prev_time);
        init_state_tracking(sys, max_E_z * unit_conversion::MHz_D_per_V_cm / calc_E_z);
        fs::create_directories(track_dir_path);
    } else {
        prev_time = log_time_at_point("Not Initializing State Tracking (force-disabled)", start_time, prev_time);
        disable_tracking();
    }

    // Stark loop
    prev_time = log_time_at_point("About to start stark loop", start_time, prev_time);

    std::vector<double> max_devs_vec(nLowestStates, -std::numeric_limits<double>::infinity());
    std::vector<int> max_devdx_vec(nLowestStates, -1);

    const double scale_Ez = max_E_z - min_E_z;
    const double scale_Ez_mhz = scale_Ez * unit_conversion::MHz_D_per_V_cm;
    const double offset_Ez_mhz = min_E_z;

    // note: we now need to go backwards to track from calc_E_z
    for (int fdx = nStarkIterations - 1; fdx >= 0; fdx--) {
        double field_divisor = nStarkIterations - 1.0;
        double Ez_fdx_mhz = (scale_Ez_mhz) * (fdx / field_divisor) + offset_Ez_mhz;
        double Ez_V_cm = Ez_fdx_mhz / unit_conversion::MHz_D_per_V_cm;
        
        // recalaculate H_tot -- from scratch to avoid accumulation of error
        // calc.H_tot.setZero();
        sys.H_tot = sys.H_rot.toDenseMatrix() + /**/ sys.H_hfs + /**/ dcomplex(Ez_fdx_mhz / calc_E_z) * sys.H_stk;

#ifdef USE_DEVONSHIRE
        sys.H_tot += sys.H_dev;
#endif
        sys.diagonalize();
        // Update tracking and then output new tracking info
        if (do_tracking) {
            update_tracking(sys, vals);
            std::string csvbas;
            if (round_files) {
                csvbas = fmt::format("{}.csv", std::lround(Ez_V_cm));
            } else {
                csvbas = fmt::format("{:.1f}.csv", Ez_V_cm);
            }
            auto track_csv_path = track_dir_path / csvbas;
            std::ofstream os(track_csv_path);
            output_tracking_info(os);
        }
        // energy output
        oEs << fmt::format("{},{}", fdx, Ez_V_cm);
        for (size_t idx = 0; idx < sys.nBasisElts; idx++) {
            oEs << fmt::format(",{}", std::real(sys.Es[idx]));
        }
        oEs << std::endl;

        // f = 0 singlet
        int32_t gnd_idx = closest_state(sys, 0);
        int32_t _if00 = gnd_idx;
        double E = std::real(sys.Es[gnd_idx]);

        // energy differences for f = 1 triplet
        std::vector<double> dEs(nLowestStates);
        std::vector<int> idxs(nLowestStates);

        // measure deviation of m_f for each of the "lowest" states
        for (int sdx = 0; sdx < nLowestStates; sdx++) {
            aef::universal_diatomic_basis_vec v = lowest_states[sdx];
            idxs[sdx] = closest_state(sys, lowest_idxs[sdx], _if00);
            double dev_mf = std::abs(expectation_values(sys, idxs[sdx]).m_f - v.m_f);

            if (dev_mf > max_devs_vec[sdx]) {
                max_devs_vec[sdx] = dev_mf;
                max_devdx_vec[sdx] = fdx;
            }

            dEs[sdx] = std::real(sys.Es[idxs[sdx]]) - E;
        }

        double stark_scale = std::nan("");
        {
            double mu_e = baf_constants::mu_e;
            (void)pCalc->get_parameter("mu_E", mu_e);
            double stark_scale = Ez_V_cm * mu_e * unit_conversion::MHz_D_per_V_cm;
        }

        std::cout << fmt::format("Electric field strength is {} V/cm, stark scale is {} MHz", Ez_V_cm, stark_scale) << std::endl;
        std::cout << fmt::format("Gnd state expectation values: {}", expectation_values(sys, gnd_idx)) << std::endl;

        for (int sdx = 0; sdx < nLowestStates; sdx++) {
            std::cout << fmt::format("#{} state expectation values: {}", sdx, expectation_values(sys, idxs[sdx])) << std::endl;
        }

        std::cout << fmt::format("Closest Energy-estate to 0-E-field gnd state is "
            "{}, with energy {}", gnd_idx, E) << std::endl;
        // write energy differences to stark log and to standard out
        {
            oStk << Ez_V_cm << "," << E;
            std::cout << stark_scale << "," << E;
            // need to start at 1 to because the lowest state is accounted for by E
            for (int sdx = 1; sdx < nLowestStates; sdx++) {
                oStk << "," << dEs[sdx];
                std::cout << "," << dEs[sdx];
            }
            oStk << std::endl;
            std::cout << std::endl;
        }

        // collect energy of the lowest group of states
        for (int sdx = 0; sdx < nLowestStates; sdx++) {
            lowest_energies(fdx, sdx) = EVAL(sys.Es[idxs[sdx]]);
        }

        std::string dev_out_fname;
        if (round_files) {
            dev_out_fname = fmt::format("info_Ez_{}.csv", std::lround(Ez_V_cm));
        } else {
            dev_out_fname = fmt::format("info_Ez_{:.1f}.csv", Ez_V_cm);
        }
        std::ofstream dout(devpath / dev_out_fname);
        output_state_info(dout, sys
#ifndef DONT_USE_CUDA
        , vals
#endif
        );
    }

#if 1
    { // print header line
        const char* sep = "";
        for (int sdx = 0; sdx < nLowestStates; sdx++) {
            std::cout << fmt::format("{}E{}", sep, sdx);
            sep = ", ";
        }
        std::cout << std::endl;
    }
    for (int fdx = 0; fdx < 101; fdx++) {
        const char* sep = "";
        for (int sdx = 0; sdx < nLowestStates; sdx++) {
            std::cout << fmt::format("{}{}", sep, lowest_energies(fdx, sdx));
            sep = ", ";
        }
    }
#endif

    std::cout << "--------- stark loop completed ---------" << std::endl;
    prev_time = log_time_at_point("Completed stark loop", start_time, prev_time);
    std::cout << fmt::format("Explicit m_f degeneracy breaking coeff is {:.4} Hz",
        aef::raf_constants::e_mf_break * 1E6) << std::endl;
    for (int sdx = 0; sdx < nLowestStates; sdx++) {
        std::cout << fmt::format("Maximum m_f deviation for {} is {} at index {}",
            lowest_states[sdx], max_devs_vec[sdx], max_devdx_vec[sdx]) << std::endl;
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
//   5. Go to Project > Add New Item to create new code files, or Project > Add
//   Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project
//   and select the .sln file

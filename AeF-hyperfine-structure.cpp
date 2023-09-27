// AeF-hyperfine-structure.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
//
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


#include "AeF-hyperfine-structure.inl"
// int32_t closest_approx(Eigen)

#undef MATRIX_ELT_DEBUG

/// <summary>
/// Calculates the expectation values
/// </summary>
/// <param name="calc">HyperfineCalculator: contains operator matrix elements
/// and states</param> <param name="E_idx">the index of Energy level to
/// calculate with</param> <returns></returns>
j_basis_vec expectation_values(HyperfineCalculator& calc, int32_t E_idx) {
    //
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    j_basis_vec out;
#ifdef _WIN32
    SecureZeroMemory((void*)&out, sizeof(j_basis_vec));
#else
    explicit_bzero((void*)&out, sizeof(j_basis_vec));
#endif
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
        j_basis_vec bs_ket = calc.basis[kidx];

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
        out.n += prob * bs_ket.n;
        out.j += prob * bs_ket.j;
        out.f += prob * bs_ket.f;
        out.m_f += prob * bs_ket.m_f;
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }

    out.n /= prob_tot;
    out.j /= prob_tot;
    out.f /= prob_tot;
    out.m_f /= prob_tot;

    return out;
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
    explicit_bzero((void*)&out, sizeof(j_basis_vec));
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
        ex_parity += prob*std::pow(-1, bs_ket.n);
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }
    return ex_parity / prob_tot;
}

static inline double diff_states(j_basis_vec v1, j_basis_vec v2) {
    constexpr double cn = 1.5;
    constexpr double cj = 0.5; // 1.0;
    constexpr double cf = 3.0;
    constexpr double cm = 4.0; // 100.0;

    double dn = (v1.n - v2.n);
    double dj = (v1.j - v2.j);
    double df = (v1.f - v2.f);
    double dm = (v1.m_f - v2.m_f);
    return cn * dn * dn + cj * dj * dj + cf * df * df + cm * dm * dm;
}

/// <summary>
/// Finds the diagonlized eignestate "closest" to the basis state with index
/// "ket_idx"
/// </summary>
/// <param name="calc">System basis and operator elements</param>
/// <param name="ket_idx">basis state</param>
/// <param name="exclude_Eidx">(optional) an energy eigenstate to "exclude" from
/// being the closest state.  Intended to fix the </param> <returns></returns>
int32_t closest_state(HyperfineCalculator& calc, int32_t ket_idx,
    int32_t exclude_Eidx = -1) {
    double chisq = (double)std::numeric_limits<double>::infinity();
    j_basis_vec ket = calc.basis[ket_idx];
    int32_t closest_idx = -1;

    for (int32_t Eidx = 0; Eidx < calc.nBasisElts; Eidx++) {
        j_basis_vec expect_ket = expectation_values(calc, Eidx);
        double localx2 = diff_states(ket, expect_ket);

        if (localx2 < chisq && Eidx != exclude_Eidx) {
            chisq = localx2;
            closest_idx = Eidx;
        }
    }

    return closest_idx;
}


time_point<system_clock> log_time_at_point(
    const char* desc,
    std::chrono::time_point<std::chrono::system_clock>& start,
    std::chrono::time_point<std::chrono::system_clock>& prev) {
    using namespace std::chrono;
    time_point<system_clock> curr_time = system_clock::now();
    std::string stime = fmt::format("{0:%F}-{0:%H%M}{0:%S}", curr_time);
    using d_seconds = std::chrono::duration<double>;
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
    return curr_time;
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
    ///////////////////// constants

    /// <summary>
    /// WARNING: the wrong value was originally used in the conversion factor
    /// This was supposed to be 50 kV/cm but is actually 500 kV/cm.
    /// </summary>
    constexpr double E_z_orig = unit_conversion::MHz_D_per_V_cm * 500 * 1000;
    /// <summary>
    /// 50 kV/cm = 25170 MHz/D is the field strength used to calculate H_stk.
    /// DO NOT CHANGE THIS.
    /// Instead use the E_max argument
    /// </summary>
    constexpr double calc_E_z = unit_conversion::MHz_D_per_V_cm * 50 * 1000;

#define USE_DEVONSHIRE

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
    // 2023-07-12: change output dir to output instead of oana
    auto dpath = fs::path("output");
    std::chrono::time_point<std::chrono::system_clock> start_time =
        std::chrono::system_clock::now();
    std::string stime = fmt::format("{0:%F}-{0:%H%M}{0:%S}", start_time);
    std::chrono::time_point<std::chrono::system_clock> prev_time = start_time;
    dpath /= stime;
    fs::create_directories(dpath);

    int param_nmax = 20;
    bool enable_debug_log = false;
    bool load_from_file = false;
    std::string loadname = "";
    bool print_extras = true;
    bool output_Es = true;
    size_t nStarkIterations = 101;
    double min_E_z = 0;
    double max_E_z = calc_E_z;

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
        std::cout << "AeF Hyperfine Structure version compiled on " << __DATE__ << " "
            << __TIME__ << ", git commit " << aef_git_commit << std::endl;
        std::cout << "Git status is " << dirty << " string {" << status << "}" << std::endl;
        std::cout << "Start time is " << start_time << std::endl;
        std::cout << fmt::format("Eigen will use {} threads", Eigen::nbThreads()) << std::endl;
    }
#ifdef _OPENMP
    std::cout << "Reconfiguring openmp to use the correct number of threads (the number of physical cores)." << std::endl;
    int num_physical_cores = get_num_cores();
    omp_set_num_threads(num_physical_cores);
    Eigen::setNbThreads(num_physical_cores);
    std::cout << fmt::format("OpenMP/Eigen will use {} threads", num_physical_cores) << std::endl;
#endif


    init_rng();

    j_basis_vec v(1, .5, 0, 0);
    double E_rot = std::real(v.H_rot());
    dcomplex H_hfs = v.H_hfs(v);
    
    dcomplex H_st = v.H_st(v, E_z_orig);
    std::string str = fmt::format("{}: E_rot={} MHz, E_hfs={} MHz, E_st(500kV/cm) = {} MHz", v,
            E_rot, H_hfs, H_st);

    std::cout << str << std::endl;

    // maximum value of the n quantum number.  There will be 4*(nmax**2) states
    int nmax = param_nmax;
    HyperfineCalculator calc(nmax, calc_E_z, K);

    std::cout << fmt::format("nmax is {}, E_z is {} MHz/D, K is {} MHz ({})",
        nmax, calc_E_z, K, devstatus)
        << std::endl;
    if (load_from_file) {
        std::string logstr = fmt::format("Loading matrix elements from {}", loadname);
        prev_time = log_time_at_point(logstr.c_str(), start_time, prev_time);
        bool result = calc.load_matrix_elts(loadname);

        if (!result) {
            std::cout << "couldn't load " << loadname << std::endl;
            exit(-1);
        }
        logstr = fmt::format("Finished loading matrix elements from {}", loadname);
        prev_time = log_time_at_point(logstr.c_str(), start_time, prev_time);
    } else {
        prev_time = log_time_at_point("Starting matrix element calculations", start_time, prev_time);
        calc.calculate_matrix_elts();
        calc.diagonalize_H();
        if (nmax >= 20)
            calc.save_matrix_elts(dpath / "matrix.dat");

        prev_time = log_time_at_point("Finished matrix elt calcs", start_time, prev_time);
    }

    if (print_extras) {
        Eigen::VectorXcd Es = calc.Es;
        std::cout << "----------- Stark-Shifted -----------" << std::endl;
        std::cout << "Level, Energy (MHz)" << std::endl;
        double EPrev = 0;
        for (int i = 0; i < calc.nBasisElts; i++) {
            double dE = std::real(Es[i]) - EPrev;
            std::cout << i << ", " << std::real(Es[i]) << "MHz, DeltaE = " << dE
                << " MHz, " << expectation_values(calc, i).ket_string()
                << std::endl;
            EPrev = std::real(Es[i]);
        }
        std::cout << std::endl << std::endl;

        std::cout << "----------- NO STARK -----------" << std::endl;
        calc.H_tot -= calc.H_stk;
        calc.diagonalize_H();

        Es = calc.Es;
        EPrev = 0;
        std::cout << "Level, ket, Energy (MHz)" << std::endl;
        for (int i = 0; i < calc.nBasisElts; i++) {
            double dE = std::real(Es[i]) - EPrev;
            std::cout << i << ", " << std::real(Es[i]) << "MHz, DeltaE = " << dE
                << " MHz, " << expectation_values(calc, i).ket_string()
                << std::endl;
            EPrev = std::real(Es[i]);
            // std::cout << i << ", " << calc.basis[i] << ", " << Es[i] << std::endl;
        }
    } else {
        calc.H_tot -= calc.H_stk;
        calc.diagonalize_H();
    }

    // create output file
    auto fpath = dpath / "stark_shift_gnd.csv";
    std::ofstream oStk(fpath, std::ios::trunc | std::ios::out);

    auto epath = dpath / "stark_spectrum.csv";
    std::ofstream oEs (epath, std::ios::trunc | std::ios::out);
    oEs << "idx,E-field (V/cm)";
    for (size_t idx = 0; idx < calc.nBasisElts; idx++) {
        oEs << fmt::format(",E{}", idx);
    }
    oEs << std::endl;

    j_basis_vec gnd = j_basis_vec::from_index(0);
    j_basis_vec f00 = gnd;
    int32_t if00 = f00.index();
    // n = 0, j = 0.5, f = 1 hyperfine triplet
    j_basis_vec f1t(0, 0.5, 1, -1);
    int32_t if1t = f1t.index();
    j_basis_vec f10(0, 0.5, 1, 0);
    int32_t if10 = f10.index();
    j_basis_vec f11(0, 0.5, 1, 1);
    int32_t if11 = f11.index();

    std::cout << "if1t=" << if1t << " if10=" << if10 << " if11=" << if11
        << std::endl;
    std::cout << fmt::format("f00={}; f1t={}, f10={}, f11={}", f00, f1t, f10, f11)
        << std::endl;

    // oStk << "E-field (V/cm), Stark-shifted Energy of " << gnd.ket_string() << "
    // (MHz)";
    oStk << "E-field (V/cm), dE_gnd" << ", dE_f1t, dE_f10, dE_f11" << std::endl;
    assert(calc.H_tot.rows() == calc.H_stk.rows());

#define USE_REAL_Es
#ifdef USE_REAL_Es
    typedef double etype;
#define EVAL(val) std::real(val) 
#else
    typedef dcomplex etype;
#define EVAL(val) val
#endif
    std::vector<etype> E0s(nStarkIterations);
    std::vector<etype> E1s(nStarkIterations);
    std::vector<etype> E2s(nStarkIterations);
    std::vector<etype> E3s(nStarkIterations);
#ifndef USE_DEVONSHIRE
    // note: devonshire potential doesn't conserve m_f
    for (int idx = 0; idx < calc.nBasisElts; idx++) {
        j_basis_vec v1 = calc.basis[idx];
        for (int jdx = 0; jdx < calc.nBasisElts; jdx++) {
            j_basis_vec v2 = calc.basis[jdx];
            double prob = std::norm(calc.H_tot(idx, jdx));
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
    std::cout << "does H_tot commute with d10? " << commutes(calc.H_tot, calc.d10)
        << std::endl;
    std::cout << "does H_tot commute with d11? " << commutes(calc.H_tot, calc.d11)
        << std::endl;
    std::cout << "does H_tot commute with d1t? " << commutes(calc.H_tot, calc.d1t)
        << std::endl;
    std::cout << std::endl;

    std::cout << "Is d10  all zero " << calc.d10.isZero(1E-6) << std::endl;
    std::cout << "Is d11  all zero " << calc.d11.isZero(1E-6) << std::endl;
    std::cout << "Is d1t  all zero " << calc.d1t.isZero(1E-6) << std::endl;
    std::cout << "Is Hdev all zero " << calc.H_dev.isZero(1E-6) << std::endl;



    // Stark loop
    prev_time = log_time_at_point("About to start stark loop", start_time, prev_time);
    double max_dev_mf_f00 = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f00 = -1;
    double max_dev_mf_f10 = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f10 = -1;
    double max_dev_mf_f1t = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f1t = -1;
    double max_dev_mf_f11 = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f11 = -1;

    const double scale_Ez = max_E_z - min_E_z;
    const double scale_Ez_mhz = scale_Ez * unit_conversion::MHz_D_per_V_cm;
    const double offset_Ez_mhz = min_E_z;

    for (int fdx = 0; fdx < nStarkIterations; fdx++) {
        double field_divisor = nStarkIterations - 1.0;
        double Ez_fdx_mhz = (scale_Ez_mhz) * (fdx / field_divisor) + offset_Ez_mhz;
#ifdef MATRIX_ELT_DEBUG
        // degenerate states will probably break this
        if (fdx == 0)
            continue;
#endif
        // recalaculate H_tot -- from scratch to avoid accumulation of error
        // calc.H_tot.setZero();
        calc.H_tot = calc.H_rot.toDenseMatrix() + /**/ calc.H_hfs + /**/ dcomplex(Ez_fdx_mhz / calc_E_z) * calc.H_stk;

#ifdef USE_DEVONSHIRE
        calc.H_tot += calc.H_dev;
#endif
        calc.diagonalize_H();

        double Ez_V_cm = Ez_fdx_mhz / unit_conversion::MHz_D_per_V_cm;
        // energy output
        oEs << fmt::format("{},{}", fdx, Ez_V_cm);
        for (size_t idx = 0; idx < calc.nBasisElts; idx++) {
            oEs << fmt::format(",{}", std::real(calc.Es[idx]));
        }
        oEs << std::endl;

        // f = 0 singlet
        int32_t gnd_idx = closest_state(calc, 0);
        int32_t _if00 = gnd_idx;
        double E = std::real(calc.Es[gnd_idx]);

        // energy differences for f = 1 triplet
        int32_t _if1t = closest_state(calc, if1t, _if00);
        double dE_f1t = std::real(calc.Es[_if1t]) - E;
        int32_t _if10 = closest_state(calc, if10, _if00);
        double dE_f10 = std::real(calc.Es[_if10]) - E;
        int32_t _if11 = closest_state(calc, if11, _if00);
        double dE_f11 = std::real(calc.Es[_if11]) - E;

        // measure deviation of m_f for each n=0,f=0 and n=0,f=1 state
#define DEC_MDEV(idx) j_basis_vec jb_f##idx = f##idx
        DEC_MDEV(00);
        DEC_MDEV(10);
        DEC_MDEV(11);
        DEC_MDEV(1t);
#define MDEV(idx)                                                              \
  do {                                                                         \
    double dev_mf_##idx =                                                      \
        std::abs(expectation_values(calc, _i##idx).m_f - jb_##idx.m_f);           \
    if (std::abs(dev_mf_##idx) > max_dev_mf_##idx) {                           \
      max_dev_mf_##idx = dev_mf_##idx;                                         \
      idx_max_mf_##idx = fdx;                                                  \
    }                                                                          \
  } while (0)

        MDEV(f00);
        MDEV(f10);
        MDEV(f1t);
        MDEV(f11);
#undef MDEV

        double stark_scale =
            Ez_V_cm * hfs_constants::mu_e * unit_conversion::MHz_D_per_V_cm;

        std::cout << fmt::format(
                "Electric field strength is {} V/cm, stark scale is {} MHz",
                Ez_V_cm, stark_scale)
            << std::endl;
        std::cout << fmt::format("Gnd state expectation values: {}",
            expectation_values(calc, gnd_idx))
            << std::endl;
        std::cout << fmt::format("f1t state expectation values: {}",
            expectation_values(calc, _if1t))
            << std::endl;
        std::cout << fmt::format("f10 state expectation values: {}",
            expectation_values(calc, _if10))
            << std::endl;
        std::cout << fmt::format("f11 state expectation values: {}",
            expectation_values(calc, _if11))
            << std::endl;

        std::cout << fmt::format("Closest Energy-estate to 0-E-field gnd state is "
            "{}, with energy {}", gnd_idx, E) << std::endl;
        oStk << Ez_V_cm << "," << E << "," << dE_f1t << "," << dE_f10 << ","
            << dE_f11 << std::endl;
        std::cout << stark_scale << "," << E << "," << dE_f1t << "," << dE_f10
            << "," << dE_f11 << std::endl;

        // collect energy of lowest three states
        E0s[fdx] = EVAL(calc.Es[0]);
        E1s[fdx] = EVAL(calc.Es[1]);
        E2s[fdx] = EVAL(calc.Es[2]);
        E3s[fdx] = EVAL(calc.Es[3]);

//#ifdef USE_DEVONSHIRE
        auto dev_out_fname = fmt::format("info_Ez_{}.csv", std::lround(Ez_V_cm));
        std::ofstream dout(devpath / dev_out_fname);
        output_state_info(dout, calc);
//#endif
    }

#if 1
    std::cout << "E0, E1, E2, E3" << std::endl;
    for (int fdx = 0; fdx < 101; fdx++) {
        std::cout << E0s[fdx] << ", " << E1s[fdx] << ", " << E2s[fdx] << ", "
            << E3s[fdx] << std::endl;
    }
#endif

    std::cout << "--------- stark loop completed ---------" << std::endl;
    prev_time = log_time_at_point("Completed stark loop", start_time, prev_time);
    std::cout << fmt::format("Explicit m_f degeneracy breaking coeff is {:.4} Hz",
        hfs_constants::e_mf_break * 1E6)
        << std::endl;
    std::cout << fmt::format("Maximum m_f deviation for {} is {} at index {}",
        f00, max_dev_mf_f00, idx_max_mf_f00)
        << std::endl;
    std::cout << fmt::format("Maximum m_f deviation for {} is {} at index {}",
        f10, max_dev_mf_f10, idx_max_mf_f10)
        << std::endl;
    std::cout << fmt::format("Maximum m_f deviation for {} is {} at index {}",
        f1t, max_dev_mf_f1t, idx_max_mf_f1t)
        << std::endl;
    std::cout << fmt::format("Maximum m_f deviation for {} is {} at index {}",
        f11, max_dev_mf_f11, idx_max_mf_f11)
        << std::endl;

    // diagonalize

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

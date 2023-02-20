// AeF-hyperfine-structure.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
//

#include <aef/aef.h>
#include <format>
#include <iostream>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <numeric>
#include <aef/teestream.hpp>

using namespace std::chrono;
namespace fs = std::filesystem;

//int32_t closest_approx(Eigen)


int32_t most_like(Eigen::MatrixXcd &d, int32_t ket_idx) {
    int32_t cdx = -1;
    double max_comp = -1;

    for (int idx = 0; idx < d.cols(); idx++) {
        dcomplex ampl = d.col(idx)(ket_idx);
        double comp = abs(ampl);

        if (comp >= max_comp) {
            cdx = idx;
            max_comp = comp;
        }
    }

    return cdx;
}

double energy_of_closest(HyperfineCalculator& calc, int32_t ket_idx) {
    int32_t bidx = most_like(calc.Vs, ket_idx);
    return std::real(calc.Es[bidx]);
}

#undef MATRIX_ELT_DEBUG

/// <summary>
/// Calculates the expectation values
/// </summary>
/// <param name="calc">HyperfineCalculator: contains operator matrix elements and states</param>
/// <param name="E_idx">the index of Energy level to calculate with</param>
/// <returns></returns>
j_basis_vec expectation_values(HyperfineCalculator &calc, int32_t E_idx) {
    //
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    j_basis_vec out;
    SecureZeroMemory((void*)&out, sizeof(j_basis_vec));
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
            throw (false);
        }

        if (prob > 0 && !set_mf) {
            expect_mf = bs_ket.m_f;
            set_mf = true;
        }
#endif
        out.n   += prob * bs_ket.n;
        out.j   += prob * bs_ket.j;
        out.f   += prob * bs_ket.f;
        out.m_f += prob * bs_ket.m_f;
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }

    out.n   /= prob_tot;
    out.j   /= prob_tot;
    out.f   /= prob_tot;
    out.m_f /= prob_tot;

    return out;
}


static inline double diff_states(j_basis_vec v1, j_basis_vec v2) {
    constexpr double cn = 1.5;
    constexpr double cj = 0.5;// 1.0;
    constexpr double cf = 3.0;
    constexpr double cm = 4.0;//100.0;
    
    double dn = (v1.n - v2.n);
    double dj = (v1.j - v2.j);
    double df = (v1.f - v2.f);
    double dm = (v1.m_f - v2.m_f);
    return cn * dn * dn + cj * dj * dj + cf * df * df + cm * dm * dm;
}

/// <summary>
/// Finds the diagonlized eignestate "closest" to the basis state with index "ket_idx"
/// </summary>
/// <param name="calc">System basis and operator elements</param>
/// <param name="ket_idx">basis state</param>
/// <param name="exclude_Eidx">(optional) an energy eigenstate to "exclude" from being the closest
/// state.  Intended to fix the </param>
/// <returns></returns>
int32_t closest_state(HyperfineCalculator& calc, int32_t ket_idx, int32_t exclude_Eidx=-1) {
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

#ifdef USE_MOST_LIKE
#define closest_state(calc, idx) most_like(calc.Vs, idx)
#elif defined(USE_CONST)
#define closest_state(calc, idx) (idx)
#endif

#define CALCULATE_HAMILTONIAN

int main() {
    // output file --> automatically make output based on current datetime
    auto dpath = fs::path("oana");
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    std::string stime = std::format("{0:%F}-{0:%H%M}{0:%S}", now);
    dpath /= stime;
    fs::create_directories(dpath);

    // create info log
    std::ofstream oLog;
    tee::teebuf* pBuf, *pErrb;
    std::streambuf *orig_coutb = std::cout.rdbuf();
    std::streambuf* orig_cerrb = std::cerr.rdbuf();

    {
        auto fpath = dpath / "out.log";
        oLog = std::ofstream(fpath, std::ios::trunc | std::ios::out);

        pBuf = new tee::teebuf(oLog.rdbuf(), orig_coutb);
        std::cout.set_rdbuf(pBuf);
        
        pErrb = new tee::teebuf(oLog.rdbuf(), orig_cerrb);
        std::cout.set_rdbuf(pErrb);
    }
    std::cout << "AeF Hyperfine Structure version compiled on " << __DATE__ << " " << __TIME__ << std::endl;


    j_basis_vec v(1, .5, 0, 0);
    double E_rot = std::real(v.H_rot());
    dcomplex H_hfs = v.H_hfs(v);
    /// <summary>
    /// WARNING: the wrong value was originally used in the conversion factor
    /// This was supposed to be 50 kV/cm but is actually 500 kV/cm.
    /// </summary>
    const double E_z_orig = unit_conversion::MHz_D_per_V_cm * 500 * 1000;
    /// <summary>
    /// Correct value of E_z is usually 50 kV/cm
    /// </summary>
    const double E_z = E_z_orig / 10.0;// *20;

    dcomplex H_st = v.H_st(v, E_z_orig);
    std::string str = std::format("{}: E_rot={} MHz, E_hfs={} MHz, E_st(500kV/cm) = {} MHz",
        v, E_rot, H_hfs, H_st);

    std::cout << str << std::endl;

    // maximum value of the n quantum number.  There will be 4*(nmax**2) states
    int nmax = 10;
    HyperfineCalculator calc(nmax, E_z);

    std::cout << std::format("nmax is {}, E_z is {} MHz/D", nmax, E_z) << std::endl;

#ifndef CALCULATE_HAMILTONIAN
    std::string spath = std::format("out/matrix_{}.dat", nmax);

    bool result = calc.load_matrix_elts(spath);

    std::cout << std::format("Loading hamiltonian from {}", spath) << std::endl;
    if (!result) {
        std::cout << "couldn't load " << spath << std::endl;
    }
#else
    std::cout << "Calculating matrix elements" << std::endl;
    calc.calculate_matrix_elts();
    calc.diagonalize_H();
#endif
#define PRINT_EXTRAS
#ifdef PRINT_EXTRAS
    Eigen::VectorXcd Es = calc.Es;
    std::cout << "----------- Stark-Shifted -----------" << std::endl;
    std::cout << "Level, Energy (MHz)" << std::endl;
    double EPrev = 0;
    for (int i = 0; i < calc.nBasisElts; i++) {
        double dE = std::real(Es[i]) - EPrev;
        std::cout << i << ", " << std::real(Es[i]) << "MHz, DeltaE = " << dE << " MHz, " << expectation_values(calc, i).ket_string() << std::endl;
        EPrev = std::real(Es[i]);
    }
    std::cout << std::endl << std::endl;
    std::cout << "----------- NO STARK -----------" << std::endl;
    //*/
    calc.H_tot -= calc.H_stk;

    calc.diagonalize_H();

    Es = calc.Es;
    EPrev = 0;
    std::cout << "Level, ket, Energy (MHz)" << std::endl;
    for (int i = 0; i < calc.nBasisElts; i++) {
        double dE = std::real(Es[i]) - EPrev;
        std::cout << i << ", " << std::real(Es[i]) << "MHz, DeltaE = " << dE << " MHz, " << expectation_values(calc, i).ket_string() << std::endl;
        EPrev = std::real(Es[i]);
        //std::cout << i << ", " << calc.basis[i] << ", " << Es[i] << std::endl;
    }
#else
    calc.H_tot -= calc.H_stk;
    calc.diagonalize_H();
#endif
    
    // create output file
    auto fpath = dpath / "stark_shift_gnd.csv";

    std::ofstream oStk(fpath, std::ios::trunc | std::ios::out);
    j_basis_vec gnd = j_basis_vec::from_index(0);
    j_basis_vec f00 = gnd; int32_t if00 = f00.index();
    // n = 0, j = 0.5, f = 1 hyperfine triplet
    j_basis_vec f1t(0, 0.5, 1, -1); int32_t if1t = f1t.index();
    j_basis_vec f10(0, 0.5, 1,  0); int32_t if10 = f10.index();
    j_basis_vec f11(0, 0.5, 1,  1); int32_t if11 = f11.index();

    std::cout << "if1t=" << if1t << " if10=" << if10 << " if11=" << if11 << std::endl;
    std::cout << std::format("f00={}; f1t={}, f10={}, f11={}", f00, f1t, f10, f11) << std::endl;

    //oStk << "E-field (V/cm), Stark-shifted Energy of " << gnd.ket_string() << " (MHz)";
    oStk << "E-field (V/cm), dE_gnd";
    oStk << ", dE_f1t, dE_f10, dE_f11" << std::endl;
    assert(calc.H_tot.rows() == calc.H_stk.rows());
#define USE_REAL_Es
#ifdef USE_REAL_Es
    double E0s[101];
    double E1s[101];
    double E2s[101];
    double E3s[101];
#else
    dcomplex E0s[101];
    dcomplex E1s[101];
    dcomplex E2s[101];
    dcomplex E3s[101];
#endif

    for (int idx = 0; idx < calc.nBasisElts; idx++) {
        j_basis_vec v1 = calc.basis[idx];
        for (int jdx = 0; jdx < calc.nBasisElts; jdx++) {
            j_basis_vec v2 = calc.basis[jdx];
            double prob = std::norm(calc.H_tot(idx, jdx));
            if (v1.m_f != v2.m_f && prob > 0) {
                DebugBreak();
                std::string ostr = std::format("ERROR v1 = {}, v2 = {}, prob = {}", v1, v2, prob);
                std::cout << ostr << std::endl;
                std::cerr << ostr << std::endl;
                assert(!(v1.m_f != v2.m_f && prob > 0));
            }
        }
    }

    double max_dev_mf_f00 = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f00 = -1;
    double max_dev_mf_f10 = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f10 = -1;
    double max_dev_mf_f1t = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f1t = -1;
    double max_dev_mf_f11 = -std::numeric_limits<double>::infinity();
    int idx_max_mf_f11 = -1;

    for (int fdx = 0; fdx < 101; fdx++) {
        double Ez_fdx = (E_z) * (fdx / 100.0);
#ifdef MATRIX_ELT_DEBUG
        // degenerate states will probably break this
        if (fdx == 0) continue;
#endif
        // recalaculate H_tot -- from scratch to avoid accumulation of error
        //calc.H_tot.setZero();
        calc.H_tot = calc.H_rot.toDenseMatrix() + /**/calc.H_hfs +/**/ dcomplex(Ez_fdx / E_z) * calc.H_stk;
        calc.diagonalize_H();

        // f = 0 singlet
        int32_t gnd_idx = closest_state(calc, 0); //most_like(calc.Vs, 0);
        int32_t _if00 = gnd_idx;
        double E = std::real(calc.Es[gnd_idx]);//energy_of_closest(calc, gnd_idx);
        double Ez_V_cm = Ez_fdx / unit_conversion::MHz_D_per_V_cm;

        // energy differences for f = 1 triplet
        int32_t _if1t = closest_state(calc, if1t);//1;// most_like(calc.Vs, if1t);
        double dE_f1t = std::real(calc.Es[_if1t]) - E;// energy_of_closest(calc, if1t) - E;
        int32_t _if10 = closest_state(calc, if10, _if00);//2;// most_like(calc.Vs, if1t);
        double dE_f10 = std::real(calc.Es[_if10]) - E;// energy_of_closest(calc, if10) - E;
        int32_t _if11 = closest_state(calc, if11);//3;// most_like(calc.Vs, if1t);
        double dE_f11 = std::real(calc.Es[_if11]) - E;//energy_of_closest(calc, if11) - E;

        // measure deviation of m_f for each n=0,f=0 and n=0,f=1 state
#define MDEV(idx) do {\
        double dev_mf_##idx = std::abs(expectation_values(calc, _i##idx).m_f - ##idx.m_f);\
        if (std::abs(dev_mf_##idx) > max_dev_mf_##idx) {\
            max_dev_mf_##idx = dev_mf_##idx;\
            idx_max_mf_##idx = fdx;\
        } }while (0) 

        MDEV(f00);
        MDEV(f10);
        MDEV(f1t);
        MDEV(f11);
#undef MDEV

        double stark_scale = Ez_V_cm * hfs_constants::mu_e * unit_conversion::MHz_D_per_V_cm;
        
        std::cout << std::format("Electric field strength is {} V/cm, stark scale is {} MHz", Ez_V_cm, stark_scale) << std::endl;
        std::cout << std::format("Gnd state expectation values: {}", expectation_values(calc, gnd_idx)) << std::endl;
        std::cout << std::format("f1t state expectation values: {}", expectation_values(calc, _if1t)) << std::endl;
        std::cout << std::format("f10 state expectation values: {}", expectation_values(calc, _if10)) << std::endl;
        std::cout << std::format("f11 state expectation values: {}", expectation_values(calc, _if11)) << std::endl;

        std::cout << std::format("Closest Energy-estate to 0-E-field gnd state is {}, with energy {}", gnd_idx, E) << std::endl;
        oStk << Ez_V_cm << "," << E << "," << dE_f1t << "," << dE_f10 << "," << dE_f11 << std::endl;
        std::cout << stark_scale << "," << E << "," << dE_f1t << "," << dE_f10 << "," << dE_f11 << std::endl;

#ifdef USE_REAL_Es
        E0s[fdx] = std::real(calc.Es[0]);
        E1s[fdx] = std::real(calc.Es[1]);
        E2s[fdx] = std::real(calc.Es[2]);
        E3s[fdx] = std::real(calc.Es[3]);
#else
        E0s[fdx] = calc.Es[0];
        E1s[fdx] = calc.Es[1];
        E2s[fdx] = calc.Es[2];
        E3s[fdx] = calc.Es[3];
#endif
    }

#if 1
    std::cout << "E0, E1, E2, E3" << std::endl;
    for (int fdx = 0; fdx < 101; fdx++) {
        std::cout << E0s[fdx] << ", " << E1s[fdx] << ", " << E2s[fdx] << ", " << E3s[fdx] << std::endl;
    }
#endif

    std::cout << "--------- stark loop completed ---------" << std::endl;

    std::cout << std::format("Explicit m_f degeneracy breaking coeff is {:.4} Hz", hfs_constants::e_mf_break * 1E6) << std::endl;
    std::cout << std::format("Maximum m_f deviation for {} is {} at index {}", f00, max_dev_mf_f00, idx_max_mf_f00) << std::endl;
    std::cout << std::format("Maximum m_f deviation for {} is {} at index {}", f10, max_dev_mf_f10, idx_max_mf_f10) << std::endl;
    std::cout << std::format("Maximum m_f deviation for {} is {} at index {}", f1t, max_dev_mf_f1t, idx_max_mf_f1t) << std::endl;
    std::cout << std::format("Maximum m_f deviation for {} is {} at index {}", f11, max_dev_mf_f11, idx_max_mf_f11) << std::endl;

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

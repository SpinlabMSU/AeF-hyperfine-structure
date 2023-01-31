// AeF-hyperfine-structure.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
//

#include <aef/aef.h>
#include <format>
#include <iostream>
#include <chrono>
#include <fstream>
#include <filesystem>

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

int main() {
    std::cout << "Hello World!\n";
    j_basis_vec v(1, .5, 0, 0);
    double E_rot = std::real(v.H_rot());
    dcomplex H_hfs = v.H_hfs(v);
    /// <summary>
    /// WARNING: the wrong value was originally used in the conversion factor
    /// This was supposed to be 50 kV/cm but is actually 500 kV/cm.
    /// </summary>
    const double E_z = unit_conversion::MHz_D_per_V_cm * 500 * 1000;

    dcomplex H_st = v.H_st(v, E_z);
    std::string str = std::format("{}: E_rot={} MHz, E_hfs={} MHz, E_st(50kV/cm) = {} MHz",
        v, E_rot, H_hfs, H_st);

    std::ofstream out("log.txt", std::ios::out);
    std::cout << str << std::endl;
    out << str << std::endl;

    // todo add code using hyperfine calculator
    int nmax = 4;
    HyperfineCalculator calc;
    std::string spath = std::format("out/matrix_{}.dat", nmax);

    bool result = calc.load_matrix_elts(spath);

    if (!result) {
        std::cout << "couldn't load " << spath << std::endl;
    }

    Eigen::VectorXcd Es = calc.Es;

    std::cout << "Level, Energy (MHz)" << std::endl;
    for (int i = 0; i < calc.nBasisElts; i++) {
        std::cout << i << ", " << Es[i] << std::endl;
    }

    std::cout << "---------- NO STARK -----" << std::endl;
    calc.H_tot -= calc.H_stk;

    calc.diagonalize_H();

    Es = calc.Es;

    std::cout << "Level, Energy (MHz)" << std::endl;
    for (int i = 0; i < calc.nBasisElts; i++) {
        std::cout << i << ", " << Es[i] << std::endl;
    }

    // find
    std::ofstream oStk(fs::path("stark_shift_gnd.csv"), std::ios::trunc | std::ios::out);
    j_basis_vec gnd = j_basis_vec::from_index(0);
    oStk << "E-field (V/cm), Stark-shifted Energy of " << gnd.ket_string() << " (MHz)" << std::endl;
    for (int fdx = 0; fdx < 101; fdx++) {
        // note that E_z is divided by 10 so that the max field is 50 kV/cm, not 500 kV/cm
        double Ez_fdx = (E_z / 10.0) * (fdx / 100.0);

        // recalaculate
        calc.H_tot += Ez_fdx / E_z * calc.H_stk;
        calc.diagonalize_H();


        int32_t idx = most_like(calc.Vs, 0);
        double E = energy_of_closest(calc, 0);
        double Ez_V_cm = Ez_fdx / unit_conversion::MHz_D_per_V_cm;

        std::cout << std::format("Closest Energy-estate to 0-E-field gnd state is {}, with energy {}", idx, E) << std::endl;
        oStk << Ez_V_cm << "," << E << std::endl;
        calc.H_tot -= Ez_fdx / E_z * calc.H_stk;
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

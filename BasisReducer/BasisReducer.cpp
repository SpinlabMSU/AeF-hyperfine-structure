// NoStark_HyperfineTester.cpp : This file contains the 'main' function. Program execution begins and ends there.
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

#include <aef/aef.h>
#include <fmt.hpp>
#include <iostream>
#include <chrono>
#include <fstream>
#include <filesystem>


using namespace std::chrono;
namespace fs = std::filesystem;

//int32_t closest_approx(Eigen)


int32_t most_like(Eigen::MatrixXcd& d, int32_t ket_idx) {
    int32_t cdx = -1;
    double max_comp = -99999;

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

    {
        constexpr int nmax = 40;
        for (int idx = 0; idx < aef::jf_basis_vec::index_of_n(40); idx++) {
            jf_basis_vec v = jf_basis_vec::from_index(idx);
            assert(v.index() == idx);
            if (v.index() != idx) {
                std::cerr << "BAD INDEX " << idx << std::endl;
                throw idx;
            }
        }
    }



    j_basis_vec v(1, .5, 0, 0);
    double E_rot = std::real(v.H_rot());
    dcomplex H_hfs = v.H_hfs(v);
    /// <summary>
    /// WARNING: the wrong value was originally used in the conversion factor
    /// This was supposed to be 50 kV/cm but is actually 500 kV/cm.
    /// </summary>
    const double E_z = unit_conversion::MHz_D_per_V_cm * 500 * 1000;

    dcomplex H_st = v.H_st(v, E_z);
    std::string str = fmt::format("{}: E_rot={} MHz, E_hfs={} MHz, E_st(50kV/cm) = {} MHz",
        v, E_rot, H_hfs, H_st);

    std::ofstream out("log.txt", std::ios::out);
    std::cout << str << std::endl;
    out << str << std::endl;

    // todo add code using hyperfine calculator
    int nmax = 1;
    HyperfineCalculator calc(nmax, 0, false);
    
#if 0
    std::string spath = fmt::format("out/matrix_{}.dat", nmax);

    bool result = calc.load_matrix_elts(spath);

    if (!result) {
        std::cout << "couldn't load " << spath << std::endl;
    }
#endif
    calc.calculate_matrix_elts();
    //calc.H_tot -= calc.H_stk;
    calc.diagonalize_H();
    Eigen::VectorXcd Es = calc.Es;
    // find
    j_basis_vec gnd = j_basis_vec::from_index(0);
    std::cout << gnd.ket_string() << std::endl;

    double E = energy_of_closest(calc, 0);

    // n = 0, j = 0.5, f = 1 hyperfine triplet
    j_basis_vec f1t(0, 0.5, 1, -1); int32_t if1t = f1t.index();
    j_basis_vec f10(0, 0.5, 1, 0); int32_t if10 = f10.index();
    j_basis_vec f11(0, 0.5, 1, 1); int32_t if11 = f11.index();

    double dE_f1t = energy_of_closest(calc, if1t) - E;
    double dE_f10 = energy_of_closest(calc, if10) - E;
    double dE_f11 = energy_of_closest(calc, if11) - E;


    std::string ostr = fmt::format("Gnd state energy: {}, Shift of f=1 m_f=-1: {}, Shift of f=1 m_f=0: {},"
        "Shift of f = 1 m_f = 1: {}", E, dE_f1t, dE_f10, dE_f11);
    std::cout << ostr << std::endl;

    dcomplex E_hfs_scalar = f10.H_hfs_scalar(f10);
    dcomplex E_hfs_tensor = f10.H_hfs_tensor(f10);

    //                             njfm_f
    std::cout << "E_hfs_scalar for |0,+,+,0>: " << E_hfs_scalar << ", E_hfs_tensor for |0, +, +, 0>" << E_hfs_tensor << std::endl;

    E_hfs_scalar = f11.H_hfs_scalar(f11);
    E_hfs_tensor = f11.H_hfs_tensor(f11);

    //                             njfm_f
    std::cout << "E_hfs_scalar for |0,+,+,1>: " << E_hfs_scalar << ", E_hfs_tensor for |0, +, +, 1>" << E_hfs_tensor << std::endl;

    E_hfs_scalar = f1t.H_hfs_scalar(f1t);
    E_hfs_tensor = f1t.H_hfs_tensor(f1t);

    //                             njfm_f
    std::cout << "E_hfs_scalar for |0,+,+,t>: " << E_hfs_scalar << ", E_hfs_tensor for |0, +, +,-1>" << E_hfs_tensor << std::endl;

    //std::cout <<  E << "," << dE_f1t << "," << dE_f10 << "," << dE_f11 << std::endl;
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

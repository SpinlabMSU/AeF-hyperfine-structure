// GenerateHamiltonianFiles.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/*
    GenerateHamiltonianFiles/GenerateHamiltonianFiles.cpp -- this program
    pre-generates matrix element data for 138BaF systems with nmax=21 to 40.
    Note: I don't actually use this very often since the matrix element
    calculation step takes only a small fraction of the time that the
    stark-loop diagonalization takes.

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

int main() {
    std::cout << "Hello World!\n";
    j_basis_vec v(1, .5, 0, 0);
    double E_rot = std::real(v.H_rot());
    dcomplex H_hfs = v.H_hfs(v);

    const double E_z = unit_conversion::BAD_MHz_D_per_V_cm * 50 * 1000;

    dcomplex H_st = v.H_st(v, E_z);
    std::string str = fmt::format("{}: E_rot={} MHz, E_hfs={} MHz, E_st(50kV/cm) = {} MHz",
        v, E_rot, H_hfs, H_st);

    std::ofstream out(fs::path(u8"log.txt"), std::ios::out);
    std::cout << str << std::endl;
    out << str << std::endl;


    // make output directory
    fs::create_directory(fs::path(u8"out"));

    // generate and write hamiltonians
    std::vector<std::chrono::microseconds> times;


    for (int nmax = 21; nmax <= 40; nmax++) {
        auto ostart = high_resolution_clock::now();
        auto start = high_resolution_clock::now();
        HyperfineCalculator calc(nmax, E_z);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        std::cout << "Constructing Hyperfine calculator with nmax = " << nmax << " took " << duration << " microseconds" << std::endl;
        out << "Constructing Hyperfine calculator with nmax = " << nmax << " took " << duration << " microseconds" << std::endl;

        start = high_resolution_clock::now();
        calc.calculate_matrix_elts();
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        std::cout << "Calculating hamiltonian (" << calc.nBasisElts << " basis elements) took " << duration << " microseconds" << std::endl;
        out << "Calculating hamiltonian (" << calc.nBasisElts << " basis elements) took " << duration << " microseconds" << std::endl;

        // diagonalization
        start = high_resolution_clock::now();
        calc.diagonalize_H();
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        std::cout << "Diagonalizing hamiltonian (" << calc.nBasisElts << " basis elements) took " << duration << " microseconds" << std::endl;
        out << "Diagonalizing hamiltonian (" << calc.nBasisElts << " basis elements) took " << duration << " microseconds" << std::endl;

        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - ostart);
        std::cout << "Whole calculation took " << duration << " microseconds" << std::endl;
        out << "Whole calculation took " << duration << " microseconds" << std::endl;
        times.push_back(duration);


        std::string outnam = fmt::format("out/matrix_{}.dat", nmax);
        calc.save_matrix_elts(outnam);
    }

    HyperfineCalculator ncalc;

    fs::path inpath(u8"out/matrix_4.dat");
    bool succ = ncalc.load_matrix_elts(inpath);
    
    Eigen::VectorXcd Es = ncalc.Es;


    std::cout << "Level, Energy (MHz)" << std::endl;
    for (int i = 0; i < ncalc.nBasisElts; i++) {
        std::cout << i << ", " << Es[i] << std::endl;
    }

    if (!succ) {
        std::cout << "ERROR" << std::endl;
    } else {
        std::cout << "SUCCESS" << std::endl;
    }

    for (int idx = 0; idx < times.size(); idx++) {
        std::cout << idx << ", " << times[idx] << std::endl;
        out << idx << ", " << times[idx] << std::endl;
    }
    out.close();
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

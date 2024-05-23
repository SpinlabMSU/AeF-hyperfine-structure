/*
    operator_visualizer.cpp : This file contains the 'main' function. Program execution begins and ends there.
   
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
#define _SILENCE_ALL_CXX23_DEPRECATION_WARNINGS
#include <fmt.hpp>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <cxxopts.hpp>
#include <filesystem>

#define private public
//#define class struct
#include <aef/aef.h>

namespace fs = std::filesystem;


TTree *generate_basis_ttree(HyperfineCalculator &calc) {
    TTree *basis_tree = new TTree("t_jbasis", "J-Basis in Tree form");
    aef::spin n, j, f, m_f;

    basis_tree->Branch("n", &n, "n/D");
    basis_tree->Branch("j", &j, "j/D");
    basis_tree->Branch("f", &f, "f/D");
    basis_tree->Branch("m_f", &m_f, "m_f/D");


    for (int idx = 0; idx < calc.nBasisElts; idx++) {
        n = calc.basis[idx].n;
        j = calc.basis[idx].j;
        f = calc.basis[idx].f;
        m_f = calc.basis[idx].m_f;
        basis_tree->Fill();
    }

    basis_tree->Write();
    return basis_tree;
}


TTree *generate_matrix_tree(HyperfineCalculator &calc, const char *name, const char *title, Eigen::MatrixXcd &op, double mag_thresh=1.0e-17, bool use_rel_thres=true) {
    TTree* matrix_tree = new TTree(name, title);

    const Eigen::Index nBasisElts = calc.nBasisElts;

    double threshold = mag_thresh;

    if (use_rel_thres) {
        double min_mag = std::numeric_limits<double>::infinity();
        double max_mag = -1;
        for (Eigen::Index jdx = 0; jdx < nBasisElts; jdx++) {
            for (Eigen::Index idx = 0; idx < nBasisElts; idx++) {
                aef::dcomplex val = op(idx, jdx);
                double mag = std::norm(val);
                min_mag = std::min(min_mag, mag);
                max_mag = std::max(max_mag, mag);
            }
        }
        threshold = mag_thresh * max_mag;
    }
    
    Eigen::Index idx, jdx;
    double mag, phase;
    aef::dcomplex val;

    matrix_tree->Branch("idx", &idx, "idx/I");
    matrix_tree->Branch("jdx", &jdx, "jdx/I");
    matrix_tree->Branch("mag", &mag, "mag/D");
    //matrix_tree->Branch("val", &val, "val/D");
    matrix_tree->Branch("phase", &phase, "phase/D");

    size_t nFills = 0;
    size_t nZeros = 0;
    size_t nFails = 0;
    size_t nElts = nBasisElts * nBasisElts;
    for (jdx = 0; jdx < nBasisElts; jdx++) {
        for (idx = 0; idx < nBasisElts; idx++) {
            val = op(idx, jdx);
            mag = std::abs(val); // note: norm --> mag squared
            phase = std::arg(val);
            // only fill
            if (mag > threshold) {
                matrix_tree->Fill();
                nFills++;
            } else if (mag > 0) {
                nFails++;
            } else {
                nZeros++;
            }
        }
    }
    std::cout << fmt::format("Filled {} entries of {} ({} %), {} zeros, {} fails", nFills, nElts, (100 * nFills) / (double)nElts, nZeros, nFails) << std::endl;
    matrix_tree->Write();
    return matrix_tree;
}


int main(int argc, char **argv) {
    std::cout << "Hello World!\n";
    std::string matrix_file_in_name = "matrix.dat";
    std::string root_file_name = "matrix.root";
    fs::path path_mat_file(matrix_file_in_name);
    fs::path absPathMat = fs::absolute(path_mat_file);
    fs::path run_path = absPathMat.parent_path();
    std::string run_name = run_path.filename().string();
    bool no_stark = false;
    
    cxxopts::Options options("operator-visualizer", "Program to prepare operator visualization ROOT files");

    options.add_options()
        ("h,help", "Print usage")
        //("e,Ez", "Electric field values to consider in V/cm (can use multiple times)", cxxopts::value<std::vector<double>>()) // allow multiple Ez
        ("d,enable_debug", "Enable debug mode", cxxopts::value<bool>()->default_value("false"))
        ("print_extras", "Print extra information", cxxopts::value<bool>()->default_value("true"))
        ("l,load", "Set file to load molecular system operators from", cxxopts::value<std::string>())
        ("o,output", "Set output ROOT file.  Default is matrix.root", cxxopts::value<std::string>())
        ("no_stark", "Remove stark potential from H_tot", cxxopts::value<bool>()->default_value("false"));
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("load")) {
        matrix_file_in_name = result["load"].as<std::string>();
    }

    if (result.count("no_stark")) {
        no_stark = result["no_stark"].as<bool>();
    }

    HyperfineCalculator calc;
    calc.load_matrix_elts(matrix_file_in_name);

    if (no_stark) {
        calc.H_tot -= calc.H_stk;
        calc.E_z = 0;
    }


    TFile rfile(root_file_name.c_str(), "RECREATE");
    (void)generate_basis_ttree(calc);
    (void)generate_matrix_tree(calc, "H_tot", "total Hamiltonian", calc.H_tot);
    (void)generate_matrix_tree(calc, "U", "Energy eigenbasis in matrix form", calc.Vs);

    // save important run information
    TDirectory* dir = gDirectory->mkdir("parameters");
    dir->cd();
    auto* param_nmax = new TParameter<double>("nmax", calc.nmax);
    param_nmax->Write();
    auto* param_enableDev = new TParameter<bool>("enableDev", calc.enableDev);
    param_enableDev->Write();
    auto* param_K = new TParameter<double>("K", calc.K);
    param_K->Write();
    auto* param_E_z = new TParameter<double>("E_z", calc.E_z);
    param_E_z->Write();
    bool stark_only = fs::exists(run_path / "STARK_ONLY.txt");
    auto* param_stark_only = new TParameter<bool>("stark_only", stark_only);
    param_stark_only->Write();
    dir->WriteObject(&(run_name), "run");
    dir->Write();
    rfile.Flush();
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

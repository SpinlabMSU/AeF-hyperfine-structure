// AeF-hyperfine-structure.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
//

#include "aef.h"
#include <format>
#include <iostream>

int main() {
  std::cout << "Hello World!\n";
  j_basis_vec v(1, .5, 0, 0);
  double E_rot = std::real(v.H_rot());
  dcomplex H_hfs = v.H_hfs(v);

  const double E_z = unit_conversion::MHz_D_per_V_cm * 50 * 1000;

  dcomplex H_st = v.H_st(v, E_z);
  std::string str =
      std::format("{}: E_rot={} MHz, E_hfs={} MHz, E_st(50kV/cm) = {} MHz", v,
                  E_rot, H_hfs, H_st);
  std::cout << str << std::endl;
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

#pragma once
#ifndef _AEF_MOLECULAR_SYSTEM_H

#define _AEF_MOLECULAR_SYSTEM_H 1
#include "Eigen/Eigen"
#include "aef.h"
#include <filesystem>
#include <istream>
#include <unordered_map>
#include <vector>

/// <summary>
/// Future V2 output
/// </summary>
class MolecularSystem {
  std::vector<j_basis_vec> basis;
  std::unordered_map<std::string, Eigen::MatrixXcd *> operators;

public:
  MolecularSystem();
  ~MolecularSystem();

  bool load(std::istream &in);
  bool save(std::ostream &out);
};

#endif //_AEF_MOLECULAR_SYSTEM_H
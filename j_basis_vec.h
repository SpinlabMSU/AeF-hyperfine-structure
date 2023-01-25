#ifndef _J_BASIS_VEC_H
#define _J_BASIS_VEC_H
#pragma once

#include "aef_types.h"
#include <format>
#include <iostream>

class j_basis_vec {
public:
  spin n;
  spin j;
  spin f;
  spin m_f;

  j_basis_vec(spin n_ = 0, spin j_ = .5, spin f_ = 0, spin m_f_ = 0)
      : n(n_), j(j_), f(f_), m_f(m_f_) {}

  int index();

  inline spin Nmag() { return n * (n + 1); }
  inline spin Jmag() { return j * (j + 1); }
  inline spin Fmag() { return f * (f + 1); }
  constexpr inline spin Imag() { return half * (half + 1); }
  constexpr inline spin Smag() { return Imag(); }

  inline spin n_dot_s() { return half * (Jmag() - Nmag() - Smag()); }

  dcomplex H_rot();
  dcomplex H_hfs(j_basis_vec other);
  dcomplex H_st(j_basis_vec other, double E_z = 1.0);

  std::string ket_string();

  static j_basis_vec from_index(int idx);

private:
  dcomplex H_hfs_scalar(j_basis_vec other);
  dcomplex H_hfs_tensor(j_basis_vec other);
};

std::ostream &operator<<(std::ostream &os, j_basis_vec &v);


template <> struct std::formatter<j_basis_vec> : std::formatter<std::string> {
  auto format(j_basis_vec v, format_context &ctx) {
    return formatter<string>::format(std::format("|n={},j={},f={},m_f={}>", v.n, v.j, v.f, v.m_f), ctx);
  }
};

#endif
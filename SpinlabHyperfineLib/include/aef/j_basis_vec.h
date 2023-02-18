#ifndef _J_BASIS_VEC_H
#define _J_BASIS_VEC_H
#pragma once

#include "aef_types.h"
#include <format>
#include <iostream>

/// <summary>
/// Represents an element of the |(ns)jifm_f> basis of ^{138}BaF's ground electronic and vibrational
/// band.  Note that in this band, s = i = 1/2, so these are somtimes omitted from the kets, which are
/// then written as |njfm_f>.
/// </summary>
class j_basis_vec {
public:
    /// <summary>
    /// Quantum number describing the "classical"/"rotational" angular momentum of the molecule.
    /// This is the sum of the "rigid-body" angular momentum of the molecule and the electron orbital
    /// angular momentum.
    /// </summary>
    spin n;
    /// <summary>
    /// Quantum number related to the coupled "rotational" and electron-spin angular momenta
    /// \vec{J} = \vec{N} + \vec{S}
    /// </summary>
    spin j;
    /// <summary>
    /// The total angular momentum of the molecule.
    /// \vec{F} = \vec{J} + \vec{I}
    /// </summary>
    spin f;
    /// <summary>
    /// The magnetic quantum number corresponding to f.
    /// </summary>
    spin m_f;

    j_basis_vec(spin n_ = 0, spin j_ = .5, spin f_ = 0, spin m_f_ = 0)
        : n(n_), j(j_), f(f_), m_f(m_f_) {}

    /// <summary>
    /// This defines a well-ordering on the basis
    /// </summary>
    /// <returns>the index corresponding to this basis element</returns>
    int index();

    /// <summary>
    /// Calculates the eigenvalue of \abs{\vec{N}}^2
    /// </summary>
    /// <returns>The squared-magnitude of \vec{N}</returns>
    inline spin Nmag() { return n * (n + 1); }
    /// <summary>
    /// Calculates the eigenvalue of \abs{\vec{J}}^2
    /// </summary>
    /// <returns>The squared-magnitude of \vec{J}</returns>
    inline spin Jmag() { return j * (j + 1); }
    /// <summary>
    /// Calculates the eigenvalue of \abs{\vec{F}}^2
    /// </summary>
    /// <returns>The squared-magnitude of \vec{F}</returns>
    inline spin Fmag() { return f * (f + 1); }
    /// <summary>
    /// Calculates the eigenvalue of \abs{\vec{I}}^2
    /// </summary>
    /// <returns>The squared-magnitude of \vec{I}</returns>
    constexpr inline spin Imag() { return half * (half + 1); }
    /// <summary>
    /// Calculates the eigenvalue of \abs{\vec{S}}^2
    /// </summary>
    /// <returns>The squared-magnitude of \vec{S}</returns>
    constexpr inline spin Smag() { return Imag(); }

    /// <summary>
    /// Calculates the eigenvalue of \vec{N}\cdot\vec{S}
    /// </summary>
    /// <returns>\vec{N}\cdot\vec{S}</returns>
    inline spin n_dot_s() { return half * (Jmag() - Nmag() - Smag()); }

    /// <summary>
    /// Computes the rotational energy of this basis ket (since H_rot is diagonal in the 
    /// |njfm_f> basis).
    /// 
    /// This also contains the explicit symmetry breaking operations
    /// </summary>
    /// <returns>The rotaional energy of this state in MHz</returns>
    dcomplex H_rot();
    /// <summary>
    /// Computes the hyperfine Hamiltonian matrix element between two |njfm_f> basis kets (since
    /// H_hfs is not diagonal in the
    /// </summary>
    /// <param name="other">The other state to compute the matrix element with</param>
    /// <returns>The hyperfine matrix element &lt;other|H_hfs|this&gt; in MHz</returns>
    dcomplex H_hfs(j_basis_vec other);
    /// <summary>
    /// 
    /// </summary>
    /// <param name="other"></param>
    /// <param name="E_z">Electric field strength in MHz/Debye</param>
    /// <returns>The stark shift matrix element &lt;other|H_st|this&gt; in MHz</returns>
    dcomplex H_st(j_basis_vec other, double E_z = 1.0);

    /// <summary>
    /// Descibes this state as a ket
    /// </summary>
    /// <returns>A string</returns>
    std::string ket_string();
    /// <summary>
    /// Returns a string that can be used in a CSV field without quoting
    /// </summary>
    /// <returns>description suitable for CSV</returns>
    std::string ket_csv_str();

    /// <summary>
    /// This is the inverse function of j_basis_vec::index.
    /// </summary>
    /// <param name="idx">The basis index</param>
    /// <returns>The basis element corresponding to the supplied index</returns>
    static j_basis_vec from_index(int idx);

    /// <summary>
    /// Returns the lowest index corresponding to a state in the given rotational band
    /// </summary>
    /// <param name="n">The rotational quantum number</param>
    /// <returns>The smallest index corresponding to a state in the given rotational band</returns>
    static int index_of_n(spin n);

public:
    /// <summary>
    /// Evaluates the scalar part of the hyperfine matrix element
    /// </summary>
    /// <param name="other"></param>
    /// <returns></returns>
    dcomplex H_hfs_scalar(j_basis_vec other);
    /// <summary>
    /// Evaluates the rank-2 tensor part of the hyperfine matrix element.
    /// </summary>
    /// <param name="other"></param>
    /// <returns></returns>
    dcomplex H_hfs_tensor(j_basis_vec other);
};

std::ostream& operator<<(std::ostream& os, j_basis_vec& v);
bool operator == (const j_basis_vec& v1, const j_basis_vec& v2);

template <> struct std::formatter<j_basis_vec> : std::formatter<std::string> {
    auto format(j_basis_vec v, format_context& ctx) {
        return formatter<string>::format(std::format("|n={},j={},f={},m_f={}>", v.n, v.j, v.f, v.m_f), ctx);
    }
};

#endif
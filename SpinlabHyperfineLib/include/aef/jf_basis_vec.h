/*
    aef/jf_basis_vec.h -- implements the jf_basis_vec class, which describes a
    basis ket of the "J-basis" |(i1(i2(sn)j)f1)fm_f> basis for 225RaF's
    rotational and hyperfine state.

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
#ifndef _AEF_JF_BASIS_VEC_H
#define _AEF_JF_BASIS_VEC_H 1

#include "aef/aef_types.h"

#pragma once
namespace aef {
    /// <summary>
    /// Represents an element of the |(ns)jif_1fm_f> basis of ^{225}RaF's ground electronic and vibrational
    /// band.  Note that in this band, s = i_1 = i_2 = 1/2, so these are somtimes omitted from the kets, which are
    /// then written as |njfm_f>.
    /// 
    /// </summary>
    class jf_basis_vec {
        constexpr static spin I_1 = half;
        constexpr static spin I_2 = half;
        constexpr static spin S = half;
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
        /// Quantum number describing the coupled j-vector and the $^{19}$F nuclear spin
        /// \vec{F_1} = \vec{J} + \vec{I_1}
        /// </summary>
        spin f1;

        /// <summary>
        /// The total angular momentum of the molecule.
        /// \vec{F} = \vec{F_1} + \vec{I_2}
        /// </summary>
        spin f;
        /// <summary>
        /// The magnetic quantum number corresponding to f.
        /// </summary>
        spin m_f;

        jf_basis_vec(spin n_ = 0, spin j_ = .5, spin f_1_ = 0, spin f_ = 0, spin m_f_ = 0)
            : n(n_), j(j_), f1(f_1_), f(f_), m_f(m_f_) {
        }

        /// <summary>
        /// This defines a well-ordering on the basis
        /// </summary>
        /// <returns>the index corresponding to this basis element</returns>
        int index();

        /// <summary>
        /// Calculates the eigenvalue of \abs{\vec{N}}^2
        /// </summary>
        /// <returns>The squared-magnitude of \vec{N}</returns>
        inline spin Nmag() {
            return n * (n + 1);
        }
        /// <summary>
        /// Calculates the eigenvalue of \abs{\vec{J}}^2
        /// </summary>
        /// <returns>The squared-magnitude of \vec{J}</returns>
        inline spin Jmag() {
            return j * (j + 1);
        }

        /// <summary>
        /// Calculates the eigenvalue of \abs{\vec{F_1}}^2
        /// </summary>
        /// <returns>The squared-magnitude of \vec{F_1}</returns>
        inline spin F1mag() {
            return f1 * (f1 + 1);
        }

        /// <summary>
        /// Calculates the eigenvalue of \abs{\vec{F}}^2
        /// </summary>
        /// <returns>The squared-magnitude of \vec{F}</returns>
        inline spin Fmag() {
            return f * (f + 1);
        }
        /// <summary>
        /// Calculates the eigenvalue of \abs{\vec{I_1}}^2
        /// </summary>
        /// <returns>The squared-magnitude of \vec{I_1}</returns>
        constexpr inline spin I1mag() {
            return half * (half + 1);
        }
        /// <summary>
        /// Calculates the eigenvalue of \abs{\vec{I_2}}^2
        /// </summary>
        /// <returns>The squared-magnitude of \vec{I_2}</returns>
        constexpr inline spin I2mag() {
            return half * (half + 1);
        }
        /// <summary>
        /// Calculates the eigenvalue of \abs{\vec{S}}^2
        /// </summary>
        /// <returns>The squared-magnitude of \vec{S}</returns>
        constexpr inline spin Smag() {
            return I1mag();
        }

        /// <summary>
    /// Calculates the eigenvalue of \vec{N}\cdot\vec{S}
    /// </summary>
    /// <returns>\vec{N}\cdot\vec{S}</returns>
        inline spin n_dot_s() {
            return half * (Jmag() - Nmag() - Smag());
        }

        /// <summary>
        /// Computes the rotational energy of this basis ket (since H_rot is diagonal in the 
        /// |njf1fm_f> basis).
        /// 
        /// This also contains the explicit symmetry breaking operation that splits m_f even at
        /// zero field (makes sure that E-estates have definite m_f).
        /// </summary>
        /// <returns>The rotaional energy of this state in MHz</returns>
        dcomplex H_rot();
        /// <summary>
        /// Computes the hyperfine Hamiltonian matrix element between two |njfm_f> basis kets.
        /// The hyperfine Hamiltonian conserves f and m_f, and only breaks n by a small amount.
        /// </summary>
        /// <param name="other">The other state to compute the matrix element with</param>
        /// <returns>The hyperfine matrix element &lt;other|H_hfs|this&gt; in MHz</returns>
        dcomplex H_hfs(jf_basis_vec other);
        /// <summary>
        /// Computes the Stark shift matrix element between two |njfm_f> basis kets.
        /// This conserves m_f, mixing n, j, and f by an amount depending on the electric field.
        /// </summary>
        /// <param name="other">The other state</param>
        /// <param name="E_z">Electric field strength in MHz/Debye</param>
        /// <returns>The stark shift matrix element &lt;other|H_st|this&gt; in MHz</returns>
        dcomplex H_st(jf_basis_vec other, double E_z = 1.0);

        /// <summary>
        /// Computes the T^1_0 matrix element of the orientation-vector operator.
        /// The orientation vector is
        /// To calculate the spherical-tensor operator components, use the
        /// Wigner-Eckart theorem
        /// </summary>
        /// <param name="other">The other state</param>
        /// <returns>The reduced matrix element &lt;other||d^1||this&gt; </returns>
        dcomplex d10(jf_basis_vec other);

        /// <summary>
        /// Computes the T^1_1 matrix element of the orientation-vector operator.
        /// </summary>
        /// <param name="other">The other state</param>
        /// <returns>The T11 matrix element &lt;other||d^1||this&gt; </returns>
        dcomplex d11(jf_basis_vec other);

        /// <summary>
        /// Computes the T^1_{-1} matrix element of the orientation-vector operator.
        /// The orientation vector is
        /// To calculate the spherical-tensor operator components, use the
        /// Wigner-Eckart theorem
        /// </summary>
        /// <param name="other">The other state</param>
        /// <returns>The reduced matrix element &lt;other||d^1||this&gt; </returns>
        dcomplex d1t(jf_basis_vec other);

        /// <summary>
        /// Evaluates the devonshire potential matrix element between this state and other.
        /// This potential doesn't conserve any of the basis quantum numbers other than I and S
        /// (which 
        /// </summary>
        /// <param name="other">The other state</param>
        /// <param name="K">The devonshire coupling constant in MHz</param>
        /// <returns>The devonshire matrix element &lt;other|H_dev|this&gt; in MHz</returns>
        dcomplex H_dev(jf_basis_vec other, double K);

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
        /// This is the inverse function of jf_basis_vec::index.
        /// </summary>
        /// <param name="idx">The basis index</param>
        /// <returns>The basis element corresponding to the supplied index</returns>
        static jf_basis_vec from_index(int idx);

        /// <summary>
        /// Returns the lowest index corresponding to a state in the given rotational band
        /// </summary>
        /// <param name="n">The rotational quantum number</param>
        /// <returns>The smallest index corresponding to a state in the given rotational band</returns>
        static int index_of_n(spin n);

        /// <summary>
        /// The following memberr functions 
        /// </summary>
    public:

        /// <summary>
        /// Evaluates the scalar part of the hyperfine matrix element for the light nucleus (19F)
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        dcomplex H_hfs_fermi_1(jf_basis_vec other);
        /// <summary>
        /// Evaluates the rank-2 tensor part of the hyperfine matrix element for the light nucleus (19F).
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        dcomplex H_hfs_dipole_1(jf_basis_vec other);

        /// <summary>
        /// Evaluates the nuclear spin-rotation part of the hyperfine matrix element for the light nucleus (19F).
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        dcomplex H_hfs_nsr_1(jf_basis_vec other);

        /// <summary>
        /// Evaluates the scalar part of the hyperfine matrix element for the heavy nucleus (225Ra)
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        dcomplex H_hfs_fermi_2(jf_basis_vec other);
        /// <summary>
        /// Evaluates the rank-2 tensor part of the hyperfine matrix element for the heavy nucleus (225Ra).
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        dcomplex H_hfs_dipole_2(jf_basis_vec other);

        /// <summary>
        /// Evaluates the nuclear spin-rotation part of the hyperfine matrix element for the heavy nucleus (225Ra)
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        dcomplex H_hfs_nsr_2(jf_basis_vec other);

        // for internal implementation use only --> using these permits not recomputing as many wigner symbols
    private:
        dcomplex H_hfs_fermi_1(jf_basis_vec other, double w6j_jpsn, double w6j_f1jpj);
        dcomplex H_hfs_dipole_1(jf_basis_vec other, double w3j_n2np, double w6j_f1jpj, double w9j_nnpsspjjp);
        dcomplex H_hfs_nsr_1(jf_basis_vec other, double w6j_jpnsnj, double w6j_f1jpj);
        dcomplex H_hfs_fermi_2(jf_basis_vec other, double w6j_jpsn, double w6j_f1pjpjf1, double w6j_ff1pf1);
        dcomplex H_hfs_dipole_2(jf_basis_vec other, double w3j_n2np, double w6j_f1pjpjf1, double w6j_ff1pf1, double w9j_nnpsspjjp);
        dcomplex H_hfs_nsr_2(jf_basis_vec other, double w6j_jpnsnj, double w6j_f1pjpjf1, double w6j_ff1pf1);
    };
};

using aef::jf_basis_vec;
std::ostream& operator<<(std::ostream& os, jf_basis_vec& v);
bool operator == (const jf_basis_vec& v1, const jf_basis_vec& v2);

template <> struct fmt::formatter<jf_basis_vec> : fmt::formatter<std::string> {
    auto format(jf_basis_vec v, format_context& ctx) const {
        return formatter<std::string>::format(fmt::format("|n={},j={},f_1={},f={},m_f={}>", v.n, v.j, v.f1, v.f, v.m_f), ctx);
    }
};
#endif //_AEF_JF_BASIS_VEC_H
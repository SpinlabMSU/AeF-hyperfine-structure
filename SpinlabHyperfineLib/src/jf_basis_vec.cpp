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
#include "pch.h"
#include <aef/raf_constants.h>
#include <aef/aef.h>
#include "aef/jf_basis_vec.h"
namespace aef {
    namespace hfs_coeff = raf_constants;


    int jf_basis_vec::index() {
        int offset = (int)(8 * n * n);
        if (j > n) {
            offset += (int)(8 * n);
        }
        if (f1 > j) {
            offset += (int)(4 * j); // guess
        }

        if (f > f1) {
            offset += (int)(2 * f1);
        }

        offset += (int)(m_f + f);

        return offset;
    }

    dcomplex jf_basis_vec::H_rot() {
        const double nsq = Nmag();
        const double jsq = Jmag();
        const double nds = n_dot_s();

        using namespace hfs_coeff;
        return B * nsq - D * nsq * nsq + (gamma + delta * nsq) * nds;
    }

    dcomplex jf_basis_vec::H_hfs(jf_basis_vec other) {
        if (f != other.f || other.m_f != m_f) return 0;

        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;

        // pre-compute wigner symbols
        const double w6j_jpsnsj = w6j(jp, S, n, S, j, 1);
        const double w6j_f1jpj = w6j(f1, jp, I_1, 1, I_1, j);
        const double w6j_f1pjpjf1 = w6j(f1p, jp, I_1, j, f1, j);
        const double w6j_ff1pf1 = w6j(f, f1p, I_2, 1, I_2, f1);
        const double w6j_jpnsnj = w6j(jp, n, S, n, j, 1); // for NSR

        const double w3j_n2np = w3j(n, 2, np, 0, 0, 0);
        const double w9j_nnpsspjjp = w9j(n, np, 2, S, S, 1, j, jp, 1);

        // H_hfs = hfs_light_nucleus + hfs_heavy_nucleus + hfs_I1dotI2

        // hyperfine structure resulting from the light nucleus $^{19}F$
        dcomplex hfs_F = this->H_hfs_fermi_1(other, w6j_jpsnsj, w6j_f1jpj) +
            this->H_hfs_dipole_1(other, w3j_n2np, w6j_f1jpj, w9j_nnpsspjjp) +
            this->H_hfs_nsr_1(other, w6j_jpsnsj, w6j_f1jpj);
        // hyperfine structure resulting from the heavy nucleus $^{225}Ra$
        dcomplex hfs_Ra = this->H_hfs_fermi_2(other, w6j_jpsnsj, w6j_f1pjpjf1, w6j_ff1pf1) +
            this->H_hfs_dipole_2(other, w3j_n2np, w6j_f1pjpjf1, w6j_ff1pf1, w9j_nnpsspjjp) +
            this->H_hfs_nsr_2(other, w6j_jpnsnj, w6j_f1pjpjf1, w6j_ff1pf1);
        // note: all sources I have seen ignore the nuclear spin-nuclear spin term
        // it's probably too small to matter under this approximation anyways
        dcomplex hfs_I1dotI2 = 0;
        return hfs_F + hfs_Ra + hfs_I1dotI2;
    }

    dcomplex jf_basis_vec::H_st(jf_basis_vec other, double E_z) {
        const spin np = other.n;
        const spin jp = other.j;
        const spin f1p = other.f1;
        const spin fp = other.f;
        const spin m_fp = other.m_f;
        using namespace hfs_coeff;
        if (m_f != m_fp) {
            return 0;
        }
        // stark scale in MHz
        const double scale = hfs_coeff::mu_e * E_z;

        /// WARNING: f1 factors are currently complete guesses
        dcomplex xi_factors = xi(f, fp) * xi(j, jp) * xi(n, np) * xi(f1, f1p);
        dcomplex threej_factors = w3j(f, 1, fp, -m_f, 0, m_f) * w3j(n, 1, np, 0, 0, 0);
        dcomplex sixj_factors = w6j(f, 1, fp, f1p, half, f1) * w6j(f1, 1, f1p, jp, half, j) * w6j(j, 1, jp, np, half, n);
        dcomplex phase = parity(1 - m_f);
        dcomplex angular = xi_factors * threej_factors * sixj_factors * phase;
        return scale * angular;
    }

    dcomplex jf_basis_vec::H_dev(jf_basis_vec other, double K) {
        const spin np = other.n;
        const spin jp = other.j;
        const spin f1p = other.f1;
        const spin fp = other.f;
        const spin m_fp = other.m_f;

        // Prefactor
        // TODO check if this is correct
        dcomplex prf = K * xi(f, fp) * xi(j, jp) * xi(n, np) * xi(f1, f1p) * parity(-m_f);
        constexpr double sqrt_5_14 = constexpr_sqrt(5.0 / 14.0);
        // 3j factors involving f and m_f
        dcomplex f3f = sqrt_5_14 * w3j(f, 4, fp, -m_f, 4, m_fp) + sqrt_5_14 * w3j(f, 4, fp, -m_f, -4, m_fp) + w3j(f, 4, fp, -m_f, 0, m_fp);
        // 3j factor involving n --> this is a constraint on n and n'
        dcomplex f3n = w3j(n, 4, np, 0, 0, 0);
        // there need to be factors involving f1 and f1p
        //dcomplex f6j = w6j(f, 4, fp, jp, half, j) * w6j(j, 4, jp, np, half, n);
        // this is a complete guess TODO WARNING
        dcomplex f6j = w6j(f, 4, fp, f1p, half, f1) * w6j(f1, 4, f1p, jp, half, j) * w6j(j, 4, jp, np, half, n);
        // TODO
        return prf * f3f * f3n * f6j;
    }

    // MDA operators
    dcomplex jf_basis_vec::d10(jf_basis_vec other) {
        // WARNING: update this from H_St once the form is corrected
        const spin np = other.n;
        const spin jp = other.j;
        const spin f1p = other.f1;
        const spin fp = other.f;
        const spin m_fp = other.m_f;
        using namespace hfs_coeff;
        /*
        if (m_f != m_fp) {
            return 0;
        }
        */
        /// WARNING: f1 factors are currently complete guesses
        dcomplex xi_factors = xi(f, fp) * xi(j, jp) * xi(n, np) * xi(f1, f1p);
        // The difference between d1q is right here --|
        // Change that factor to q for d1q            v
        dcomplex threej_factors = w3j(f, 1, fp, -m_f, 0, m_f) * w3j(n, 1, np, 0, 0, 0);
        dcomplex sixj_factors = w6j(f, 1, fp, f1p, half, f1) * w6j(f1, 1, f1p, jp, half, j) * w6j(j, 1, jp, np, half, n);
        dcomplex phase = parity(1 - m_f);
        return xi_factors * threej_factors * sixj_factors * phase;
    }

    dcomplex jf_basis_vec::d11(jf_basis_vec other) {
        // WARNING: update this from H_St once the form is corrected
        const spin np = other.n;
        const spin jp = other.j;
        const spin f1p = other.f1;
        const spin fp = other.f;
        const spin m_fp = other.m_f;
        using namespace hfs_coeff;

        /// WARNING: f1 factors are currently complete guesses
        dcomplex xi_factors = xi(f, fp) * xi(j, jp) * xi(n, np) * xi(f1, f1p);
        // The difference between d1q is right here --|
        // Change that factor to q for d1q            v
        dcomplex threej_factors = w3j(f, 1, fp, -m_f, 1, m_f) * w3j(n, 1, np, 0, 0, 0);
        dcomplex sixj_factors = w6j(f, 1, fp, f1p, half, f1) * w6j(f1, 1, f1p, jp, half, j) * w6j(j, 1, jp, np, half, n);
        dcomplex phase = parity(1 - m_f);
        return xi_factors * threej_factors * sixj_factors * phase;
    }

    dcomplex jf_basis_vec::d1t(jf_basis_vec other) {
        // WARNING: update this from H_St once the form is corrected
        const spin np = other.n;
        const spin jp = other.j;
        const spin f1p = other.f1;
        const spin fp = other.f;
        const spin m_fp = other.m_f;
        using namespace hfs_coeff;

        /// WARNING: f1 factors are currently complete guesses
        dcomplex xi_factors = xi(f, fp) * xi(j, jp) * xi(n, np) * xi(f1, f1p);
        // The difference between d1q is right here --|
        // Change that factor to q for d1q            v
        dcomplex threej_factors = w3j(f, 1, fp, -m_f, -1, m_f) * w3j(n, 1, np, 0, 0, 0);
        dcomplex sixj_factors = w6j(f, 1, fp, f1p, half, f1) * w6j(f1, 1, f1p, jp, half, j) * w6j(j, 1, jp, np, half, n);
        dcomplex phase = parity(1 - m_f);
        return xi_factors * threej_factors * sixj_factors * phase;
    }


    std::string jf_basis_vec::ket_string() {
        return fmt::format("|n={}, j={}, f1={}, f={}, m_f={}>", n, j, f1, f, m_f);
    }

    std::string jf_basis_vec::ket_csv_str() {
        return fmt::format("|n={} j={} f1={} f={} m_f={}>", n, j, f1, f, m_f);
    }

    jf_basis_vec jf_basis_vec::from_index(int idx) {
        int n = (int)sqrt(idx / 8.0);
        idx -= 8 * n * n;

        spin j;
        if (idx >= 8 * n) {
            j = n + half;
            idx -= 8 * n;
        } else {
            j = n - half;
        }

        spin f1;
        if (idx >= 4 * j) {
            f1 = j + 0.5;
            idx -= (int)(2 * j);
        } else {
            f1 = j - 0.5;
        }
        spin f = 0;

        if (idx >= 2 * f1) {
            f = f1 + half;
        } else {
            idx = f1 - half;
        }

        spin m_f = idx - f;

        return jf_basis_vec(n, j, f1, f, m_f);
    }

    int jf_basis_vec::index_of_n(spin n) {
        return (int)(8 * n * n);
    }

    // HFS parts
    // Formulas here taken from J. Chem Phys 71, 389 (1982) [https://doi.org/10.1016/0301-0104(82)85045-3]
    dcomplex jf_basis_vec::H_hfs_fermi_1(jf_basis_vec other) {
        // Formulas here taken from J. Chem Phys 71, 389 (1982) [https://doi.org/10.1016/0301-0104(82)85045-3]
        if (n != other.n || f1 != other.f1 || f != other.f) return 0;
        constexpr double b_Fermi_F = hfs_coeff::b_F + hfs_coeff::c_F / 3.0;
        const double jp = other.j;
        spin t = n + S + j + jp + I_1 + f1 + 1;
        constexpr double mag_coeff = 3.0 / 2.0; // = q_mag(S)* q_mag(I_1); // switched to avoid numerical error
        double coeff = b_Fermi_F * mag_coeff * xi_prime(j, jp);
        double angular = w6j(jp, S, n, S, j, 1) * w6j(f1, jp, I_1, 1, I_1, j);
        return parity(t) * coeff * angular;
    }

    dcomplex jf_basis_vec::H_hfs_fermi_1(jf_basis_vec other, double w6j_jpsn, double w6j_f1jpj) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;
        using namespace hfs_coeff;

        if (n != np || f1 != f1p || f != fp) return 0.0;
        constexpr double b_Fermi_F = b_F + c_F / 3.0;
        spin t = n + S + j + jp + I_1 + f1 + 1;
        constexpr double mag_coeff = b_Fermi_F * q_mag(S) * q_mag(I_1);
        double coeff = mag_coeff * xi_prime(j, jp);
        return parity(t) * coeff * w6j_jpsn * w6j_f1jpj;
    }

    dcomplex jf_basis_vec::H_hfs_dipole_1(jf_basis_vec other) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;

        if (f != fp || f1 != f1p) return 0.0;
        double w3j_n2np = w3j(n, 2, np, 0, 0, 0);
        double w6j_f1jpj = w6j(f1, jp, I_1, 1, I_1, j);
        double w9j_nnpsspjjp = w9j(n, np, 2, S, S, 1, j, jp, 1);

        return this->H_hfs_dipole_1(other, w3j_n2np, w6j_f1jpj, w9j_nnpsspjjp);
    }

    dcomplex jf_basis_vec::H_hfs_dipole_1(jf_basis_vec other, double w3j_n2np, double w6j_f1jpj, double w9j_nnpsspjjp) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;
        using namespace hfs_coeff;

        if (f1 != f1p || f != fp) return 0.0;
        constexpr double coeff = c_F * constexpr_sqrt(10.0 / 3.0) * q_mag(S) * q_mag(I_1);
        double mag = coeff * xi_prime(n, np) * xi_prime(j, jp);
        const spin t = n + jp + I_1 + f1 + 1;
        return mag * parity(t) * w3j_n2np * w6j_f1jpj * w9j_nnpsspjjp;
    }

    dcomplex jf_basis_vec::H_hfs_nsr_1(jf_basis_vec other) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;

        if (n != np || f1 != f1p || f != fp) return 0.0;
        const double w6j_jpnsnj = w6j(jp, n, S, n, j, 1);
        const double w6j_f1jpj = w6j(f1, jp, I_1, 1, I_1, j);
        return this->H_hfs_nsr_1(other, w6j_jpnsnj, w6j_f1jpj);
    }

    dcomplex jf_basis_vec::H_hfs_nsr_1(jf_basis_vec other, double w6j_jpnsnj, double w6j_f1jpj) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;
        using namespace hfs_coeff;

        if (n != np || f1 != f1p || f != fp) return 0.0;
        const spin t = n + S + 2 * jp + I_1 + f1 + 1;
        const double mag = constexpr_sqrt(10.0 / 3.0) * c_I_F * xi_prime(j, jp) * q_mag(n) * q_mag(I_1);

        return parity(t) * mag * w6j_jpnsnj * w6j_f1jpj;
    }

    // heavy nucleus
    dcomplex jf_basis_vec::H_hfs_fermi_2(jf_basis_vec other) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;

        if (n != np || f != fp) return 0;

        const double w6j_jpsnsj = w6j(jp, S, n, S, j, 1);
        const double w6j_f1pjpjf1 = w6j(f1p, jp, I_1, j, f1, j);
        const double w6j_ff1pf1 = w6j(f, f1p, I_2, 1, I_2, f1);

        return this->H_hfs_fermi_2(other, w6j_jpsnsj, w6j_f1pjpjf1, w6j_ff1pf1);
    }

    dcomplex jf_basis_vec::H_hfs_fermi_2(jf_basis_vec other, double w6j_jpsn, double w6j_f1pjpjf1, double w6j_ff1pf1) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;
        using namespace hfs_coeff;

        if (n != np || f != fp) return 0;

        constexpr double b_Fermi_Ra = b_Ra + c_Ra / 3.0;
        spin t = n + S + 2 * j + I_1 + I_2 + 2 * f1p + f;
        constexpr double mag_coeff = b_Fermi_Ra * q_mag(S) * q_mag(I_2);
        double mag = mag_coeff * xi_prime(f1, f1p) * xi_prime(j, jp);

        return parity(t) * mag * w6j_jpsn * w6j_f1pjpjf1 * w6j_ff1pf1;
    }

    dcomplex jf_basis_vec::H_hfs_dipole_2(jf_basis_vec other, double w3j_n2np, double w6j_f1pjpjf1, double w6j_ff1pf1, double w9j_nnpsspjjp) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;
        using namespace hfs_coeff;

        if (f != fp) return 0.0;
        const spin t = n + j + 2 * f1p + I_1 + I_2 + f;
        constexpr double mag_coeff = c_Ra * constexpr_sqrt(10.0 / 3.0) * q_mag(S) * q_mag(I_2);
        const double mag = mag_coeff * xi_prime(n, np) * xi_prime(j, jp) * xi_prime(f1, f1p);

        return parity(t) * w3j_n2np * w6j_f1pjpjf1 * w6j_ff1pf1 * w9j_nnpsspjjp;
    }

    dcomplex jf_basis_vec::H_hfs_dipole_2(jf_basis_vec other) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;
        using namespace hfs_coeff;

        if (f != fp) return 0.0;

        const double w3j_n2np = w3j(n, 2, np, 0, 0, 0);
        const double w6j_f1pjpjf1 = w6j(f1p, jp, I_1, j, f1, j);
        const double w6j_ff1pf1 = w6j(f, f1p, I_2, 1, I_2, f1);
        const double w9j_nnpsspjjp = w9j(n, np, 2, S, S, 1, j, jp, 1);

        return this->H_hfs_dipole_2(other, w3j_n2np, w6j_f1pjpjf1, w6j_ff1pf1, w9j_nnpsspjjp);
    }

    dcomplex jf_basis_vec::H_hfs_nsr_2(jf_basis_vec other) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;

        if (n != np || f != fp) return 0.0;

        const double w6j_jpnsnj = w6j(jp, n, S, n, j, 1);
        const double w6j_f1pjpjf1 = w6j(f1p, jp, I_1, j, f1, j);
        const double w6j_ff1pf1 = w6j(f, f1p, I_2, 1, I_2, f1);

        return this->H_hfs_nsr_2(other, w6j_jpnsnj, w6j_f1pjpjf1, w6j_ff1pf1);
    }

    dcomplex jf_basis_vec::H_hfs_nsr_2(jf_basis_vec other, double w6j_jpnsnj, double w6j_f1pjpjf1, double w6j_ff1pf1) {
        const double np = other.n;
        const double jp = other.j;
        const double f1p = other.f1;
        const double fp = other.f;
        using namespace hfs_coeff;

        if (n != np || f != fp) return 0;
        const spin t = n + S + j + jp + I_1 + I_2 + 2 * f1p + f;
        const double mag = c_I_Ra * q_mag(I_2) * q_mag(n) * xi_prime(j, jp) * xi_prime(f1, f1p);

        return parity(t) * mag * w6j_jpnsnj * w6j_f1pjpjf1 * w6j_ff1pf1;
    }
};
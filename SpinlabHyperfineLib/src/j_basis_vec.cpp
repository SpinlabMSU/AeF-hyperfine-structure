#include "pch.h"
#include <format>

#define MATRIX_ELEMENT_DEBUG
#undef USE_ANALYTICAL_FORMULAS

int j_basis_vec::index() {
    int offset = (int)(4 * n * n);

    if (j > n) {
        offset += (int)(4 * n);
    }

    if (f > j) {
        offset += (int)(2 * j);
    }
    offset += (int)(m_f + f);
    return offset;
}


dcomplex j_basis_vec::H_rot() {
    const double nsq = Nmag();
    const double jsq = Jmag();
    const double nds = n_dot_s();

    const double gamma = hfs_constants::gamma;
    using namespace hfs_constants;

    // explicit symmetry breaking term:
    // lifts m_f degeneracy, making state-identification easier
    // note: this should be much smaller than the other terms -- for example,
    // with the default coefficient value of 0.1 Hz, this is smaller
    // than the uncertainties on the other terms.
    dcomplex handE = hfs_constants::e_mf_break * (m_f + f);

    return B * nsq - D * nsq * nsq + (gamma + delta * nsq) * nds + handE;
}


dcomplex j_basis_vec::H_hfs_scalar(j_basis_vec other) {
    // calculates the scalar part of the hyperfine hamiltonian coupling nuclear
    // spin to electron spin -- contains 2 contributions 

    if (n != other.n || f != other.f || m_f != other.m_f) {
        return 0;
    }
    const spin jp = other.j;
    // the effective coefficient contains a contribution both from the
    // \vec{I}\cdot\vec{S} term (b) and from the I_z*S_z term (c/3)
    const double eff_b = hfs_constants::b + hfs_constants::c / 3.0;
    const double prf = 3.0 / 4.0 * eff_b * xi_prime(j, jp);
    auto hg0 = w6j(half, half, 0, n, f, jp) * w6j(half, half, 0, n, f, j);
    auto hg1 = w6j(half, half, 1, n, f, jp) * w6j(half, half, 1, n, f, j);

    const dcomplex retval = prf * (hg1 - hg0);

#if defined(MATRIX_ELEMENT_DEBUG)
    if (std::abs(retval) > 1e-3) {
        std::cout << "Nonzero H_hfs_scalar for " << this->ket_string() << " and " <<
            other.ket_string() << " H_hfs_scalar = " << retval << std::endl;
    }
#endif

    return retval;
}

dcomplex j_basis_vec::H_hfs_tensor(j_basis_vec s2) {
    const spin jp = s2.j;
    const spin fp = s2.f;
    const spin m_fp = s2.m_f;
    const spin np = s2.n;
    if (m_f != m_fp || f != fp) return 0;

    const spin i = half, s = half;

#ifndef USE_ANALYTICAL_FORMULAS
    // Formulas here taken from J. Chem Phys 71, 389 (1982) [https://doi.org/10.1016/0301-0104(82)85045-3]
    // For some reason, using this causes states of the same f (and n/j**) with different m_f to remain
    // degenerate when an external electric field is applied.  This is probably the result of an implementation
    // error that I have not found.  In the mean time, don't use this.

    // ** technically, neither n nor j is a good quantum number, but they're approximately good
    // in the zero applied E-field limit.
    dcomplex retval = 0;

    constexpr dcomplex coeff = 3.0 / 2.0 * hfs_constants::c * constexpr_sqrt(10.0 / 3.0);
    dcomplex prf = coeff * xi_prime(n, np) * xi_prime(j, jp) * parity(1 + n + i + jp + f);
    dcomplex fact3j6j = w3j(np, 2, n, 0, 0, 0) * w6j(f, jp, half, 1, half, j);
    dcomplex fact9j = w9j(n, np, 2, half, half, 1, j, jp, 1);

    retval = prf * fact3j6j * fact9j;
#else
    // Formulas taken from J. Chem Phys 105, 7412 (1996)
    // this branch seems to approximately match the plots from PRA 98, 032513 (2018)
    dcomplex retval = 0.0;
    using hfs_constants::c;
    if (f == n + 1 && fp == np - 1 && np == n + 2) {
        retval = c / 2 * sqrt((n + 1.0) * (n + 2.0)) / (2.0 * n + 3.0);
    }

    if (f == n - 1 && fp == np + 1 && np == n - 2) {
        retval = c / 2 * sqrt((n - 1.0) * (n)) / (2.0 * n - 1.0);
    }

#endif
#if defined(MATRIX_ELEMENT_DEBUG)
    if (std::abs(retval) > 1e-3) {
        std::cout << "Nonzero H_hfs_tensor for " << this->ket_string() << " and " <<
            s2.ket_string() << " H_hfs_tensor = " << retval << std::endl;
        if (0) DebugBreak();
    }
#endif

    return retval;
}

dcomplex j_basis_vec::H_hfs(j_basis_vec other) {
    return H_hfs_scalar(other) + H_hfs_tensor(other);
}

dcomplex j_basis_vec::H_st(j_basis_vec other, double E_z) {
    if (m_f != other.m_f) {
        return 0;
    }

    const spin np = other.n;
    const spin jp = other.j;
    const spin fp = other.f;
    const spin m_fp = other.m_f;

    dcomplex xi_factors = xi(f, fp) * xi(j, jp) * xi(n, np);
    dcomplex threej_factors = w3j(f, 1, fp, -m_f, 0, m_f) * w3j(n, 1, np, 0, 0, 0);
    dcomplex sixj_factors = w6j(f, 1, fp, jp, half, j) * w6j(j, 1, jp, np, half, n);
    dcomplex phase = parity(1 - m_f);

    using hfs_constants::mu_e;
    dcomplex retval = mu_e * E_z * xi_factors * threej_factors * sixj_factors * phase;

    if (std::abs(retval) > 1) {
        std::cout << std::format("{} H_st {} nonzero {} MHz", *this, other, retval) << std::endl;
    }

    return retval;
}

dcomplex j_basis_vec::H_dev(j_basis_vec other, double K) {
    const spin np = other.n;
    const spin jp = other.j;
    const spin fp = other.f;
    const spin m_fp = other.m_f;

    dcomplex prf = K * xi(f, fp) * xi(j, jp) * xi(n, np) * parity(-m_f);
    constexpr double sqrt_5_14 = constexpr_sqrt(5.0 / 14.0);
    // 3j factors involving f and m_f
    dcomplex f3f = sqrt_5_14 * w3j(f, 4, fp, -m_f, 4, m_fp) + sqrt_5_14 * w3j(f, 4, fp, -m_f, -4, m_fp) + w3j(f, 4, fp, -m_f, 0, m_fp);
    // 3j factor involving n --> this is a constraint on n and n'
    dcomplex f3n = w3j(n, 4, np, 0, 0, 0);
    // 
    dcomplex f6j = w6j(f, 4, fp, jp, half, j) * w6j(j, 4, jp, np, half, n);
    return prf * f3f * f3n * f6j;
}

std::string j_basis_vec::ket_string() {
    return std::format("|n={}, j={}, f={}, m_f={}>", n, j, f, m_f);
}

std::string j_basis_vec::ket_csv_str() {
    return std::format("|n={} j={} f={} m_f={}>", n, j, f, m_f);
}

j_basis_vec j_basis_vec::from_index(int idx) {
    int n = (int)sqrt(idx / 4.0);
    idx -= 4 * n * n;

    spin j;
    if (idx >= 4 * n) {
        j = n + half;
        idx -= 4 * n;
    } else {
        j = n - half;
    }
    spin f = 0;

    if (idx >= 2 * j) {
        f = j + 0.5;
        idx -= (int)(2 * j);
    } else {
        f = j - 0.5;
    }

    spin m_f = idx - f;
    return j_basis_vec(n, j, f, m_f);
}

int j_basis_vec::index_of_n(spin n) {
    return (int)(4 * n * n);
}

std::ostream& operator<<(std::ostream& os, j_basis_vec& v) {
    return (os << std::format("|n={}, j={}, f={}, m_f={}>", v.n, v.j, v.f, v.m_f));
}

bool operator==(const j_basis_vec& v1, const j_basis_vec& v2) {
    if (v1.n != v2.n) return false;
    if (v1.j != v2.j) return false;
    if (v1.f != v2.f) return false;
    if (v1.m_f != v2.m_f) return false;
    return true;
}

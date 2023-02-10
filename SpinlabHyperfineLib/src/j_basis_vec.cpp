#include "pch.h"
#include <format>

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

    return B * nsq - D * nsq * nsq + (gamma + delta * nsq) * nds;
}


// See 
dcomplex j_basis_vec::H_hfs_scalar(j_basis_vec other) {
    // calculates \vec{I}\cdot\vec{S} part of hfs hamiltonian
    if (n != other.n or f != other.f or m_f != other.m_f) {
        return 0;
    }
    const spin jp = other.j;
    const double coeff = 3.0 / 4.0 * hfs_constants::b * xi_prime(j, jp);
    auto hg0 = w6j(half, half, 0, n, f, jp) * w6j(half, half, 0, n, f, j);
    auto hg1 = w6j(half, half, 1, n, f, jp) * w6j(half, half, 1, n, f, j);

    const dcomplex retval = coeff * (hg1 - hg0);

#if 1
    if (std::abs(retval) > 1e-3) {
        std::cout << "Nonzero H_hfs_scalar for " << this->ket_string() << " and " <<
            other.ket_string() << " H_hfs_scalar = " << retval << std::endl;
    }
#endif

    return retval;
}

#define USE_3J
#define NEW_HFS_TENSOR
//#define MULTI_M_IS_TEST

#ifdef USE_3J

#ifndef NEW_HFS_TENSOR
static inline dcomplex hfs_tensor_angular(double j, double m_j, double jp,
    double m_jp, double f, double m_f,
    double fp, double m_fp, double n,
    double m_n, double np, double m_np, 
    double m_i, double m_s) {
    const spin i = half, s = half;
    dcomplex w0 = w3j(i, jp, fp, m_i, m_jp, -m_fp);
    dcomplex w1 = w3j(s, np, jp, m_s, m_np, -m_jp);
    dcomplex w2 = w3j(i, j, f, m_i, m_j, -m_f);
    dcomplex w3 = w3j(s, n, j, m_s, m_n, -m_j);
    dcomplex retval = w0 * w1 * w2 * w3 * m_i * m_s * parity(-m_j - m_jp);

#if 0
    if (std::abs(retval) >= 1e-5) {
        DebugBreak();
        std::cout << "NONZERO" << std::endl;
    }
#endif

    return retval;
}
#else
static inline dcomplex hfs_tensor_angular(
    double j, double jp, double f, 
    double m_f, double fp, double m_fp,
    double n, double np, double m_i, double m_s) {
    //

    const spin m_j = m_f - m_i;
    const spin m_jp = m_fp - m_i;
    const spin m_n = m_j - m_s;
    const spin m_np = m_jp - m_s;

    const spin s = half;
    const spin i = half;

    dcomplex snj  = w3j(n,  s, j , m_n,  m_s, -m_j);
    dcomplex snjp = w3j(np, s, jp, m_np, m_s, -m_jp);

    dcomplex isf  = w3j(i, j , f , m_i, m_j , -m_f);
    dcomplex isfp = w3j(i, jp, fp, m_i, m_jp, -m_fp);

//#define MATCH_PAPER
#ifdef MATCH_PAPER
    // doing this seems to match PRA 98, 032513 (2018) --> not at nonzero stark shifts
    dcomplex par = parity(1 - f + n + np - m_j - m_jp - m_f - m_fp);

    return par * snj * isf * isfp * snjp / 6.0;
#else
    // this is the actual result I calculated
    //dcomplex par = parity(n + np - m_j - m_jp - m_f - m_fp);
    dcomplex par = parity(-n - np - m_j - m_jp - m_f - m_fp);

    return par * snj * isf * isfp * snjp;
#endif
}
#endif
#else
static inline dcomplex hfs_tensor_angular(
    double j, double jp, 
    double f, double m_f, double fp, double m_fp,
    double n, double np, double m_i, double m_s
#ifdef MULTI_M_IS_TEST
    , double m_ip, double m_sp) {
#else
    ){
    const spin m_ip = m_i;
    const spin m_sp = m_s;
#endif
    dcomplex prf = m_i * m_s;

    const spin m_j  = m_f  - m_i;
    const spin m_n  = m_j  - m_s;
    const spin m_jp = m_fp - m_ip;
    const spin m_np = m_jp - m_sp;

    const spin s = half;
    const spin i = half;

    double cg  = cg_coeff(i, m_i , j , m_j , f , m_f ) * cg_coeff(s, m_s , n , m_n , j , m_j );
    double cgp = cg_coeff(i, m_ip, jp, m_jp, fp, m_fp) * cg_coeff(s, m_sp, np, m_np, jp, m_jp);

    return m_i * m_s * cg * cgp;
}

#endif


dcomplex j_basis_vec::H_hfs_tensor(j_basis_vec s2) {
    /*
      j = self.j; jp = other.j; f = self.f; fp = other.f
      m_f = self.m_f; m_fp = other.m_f
      n = self.n

      if self.n != other.n or self.f != other.f or self.m_f != other.m_f:
        return 0

      htensor = 0
      prf = c * abs(xi(f, fp)) * xi(j, jp) * parity(-2*m_f)
      i = s = half
      for m_j in np.arange(-j, j+1, 1):
          for m_jp in np.arange(-jp, jp+1, 1):
            for m_i in (-half, half):
              for m_s in (-half, half):
                for m_n in range(-n, n+1, 1):
                  ##
                  w0 = w3j(i,jp,fp,m_i,m_jp,-m_fp)
                  w1 = w3j(s,n,jp,m_s,m_n,-m_jp)
                  w2 = w3j(i,j,f,m_i,m_j,-m_f)
                  w3 = w3j(s,n,j,m_s,m_n,-m_j)
                  htensor += w0*w1*w2*w3*m_i*m_s*parity(-m_j-m_jp)
       return prf * htensor
      */
    const spin jp = s2.j;
    const spin fp = s2.f;
    const spin m_fp = s2.m_f;
    const spin np = s2.n;
    if (m_f != m_fp || f != fp) return 0;
    if (n != s2.n) return 0;

    const spin i = half, s = half;

    dcomplex prf = hfs_constants::c;// *xi(jp, j)* xi_prime(f, fp);
#ifdef USE_3J
    prf *= xi(jp, j) * xi_prime(f, fp);
#endif

    dcomplex htensor = 0;

#if not defined(USE_3J) and defined(MULTI_M_IS_TEST)
    for (int comb = 0; comb < (1 << 4); comb++) {
        const spin m_i = (comb & (1 << 0)) ? -half : half;
        const spin m_s = (comb & (1 << 1)) ? -half : half;

        const spin m_ip = (comb & (1 << 2)) ? -half : half;
        const spin m_sp = (comb & (1 << 3)) ? -half : half;
        htensor += hfs_tensor_angular(j, jp, f, m_f, fp, m_fp, n, np, m_i, m_s, m_ip, m_sp);
    }
#elif defined(NEW_HFS_TENSOR)
    // args are:                 (j, jp, f, m_f, fp, m_fp, n, np,  m_i,   m_s )
    htensor += hfs_tensor_angular(j, jp, f, m_f, fp, m_fp, n, np, -half, -half);
    htensor += hfs_tensor_angular(j, jp, f, m_f, fp, m_fp, n, np, -half, +half);
    htensor += hfs_tensor_angular(j, jp, f, m_f, fp, m_fp, n, np, +half, -half);
    htensor += hfs_tensor_angular(j, jp, f, m_f, fp, m_fp, n, np, +half, +half);
#else
    for (int dm_j = (int)(-2 * j); dm_j <= (int)(2 * j); dm_j += 2) {
        const double m_j = dm_j / 2.0;
        for (int dm_jp = (int)(-2 * j); dm_jp <= (int)(2 * j); dm_jp += 2) {
            const double m_jp = dm_jp / 2.0;
            for (double m_n = -n; m_n <= n; m_n++) {
                for (double m_np = -np; m_np <= np; m_np++) {
                    // manually unrolled innermost loops
                    htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, fp, m_fp, n, m_n, np, m_np,
                        -half, -half);
                    htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, fp, m_fp, n, m_n, np, m_np,
                        -half, +half);
                    htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, fp, m_fp, n, m_n, np, m_np,
                        +half, -half);
                    htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, fp, m_fp, n, m_n, np, m_np,
                        +half, +half);
                }
            }
        }
    }
#endif
    dcomplex retval = prf * htensor;
#if 0
    if (std::abs(htensor) > 1e-3) {
        std::cout << "Nonzero htensor for " << this->ket_string() << " and " << s2.ket_string() << " htensor = " << htensor << std::endl;
        DebugBreak();
    }
#endif

#if 1
    if (std::abs(retval) > 1e-3) {
        std::cout << "Nonzero H_hfs_tensor for " << this->ket_string() << " and " << 
            s2.ket_string() << " H_hfs_tensor = " << retval << std::endl;
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

    dcomplex xi_factors = xi(f, other.f) * xi(j, jp) * xi(n, np) * xi(m_f, m_fp);

    dcomplex threej_factors =
        w3j(f, 1, fp, -m_f, 0, m_f) * w3j(n, 1, np, 0, 0, 0);

    dcomplex sixj_factors =
        w6j(f, 1, fp, jp, half, j) * w6j(j, 1, jp, np, half, n);

    dcomplex phase = parity(1 - m_f);

    return hfs_constants::mu_e * E_z * xi_factors * threej_factors *
        sixj_factors * phase;
}

std::string j_basis_vec::ket_string() {
    return std::format("|n={}, j={}, f={}, m_f={}>", n, j, f, m_f);
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
    return (
        os << std::format("|n={}, j={}, f={}, m_f={}>", v.n, v.j, v.f, v.m_f));
}

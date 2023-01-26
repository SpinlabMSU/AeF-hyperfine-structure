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

dcomplex j_basis_vec::H_hfs_scalar(j_basis_vec other) {
    // calculates \vec{I}\cdot\vec{S} part of hfs hamiltonian
    if (n != other.n or f != other.f or m_f != other.m_f) {
        return 0;
    }
    const spin jp = other.j;
    const double coeff = 3.0 / 4.0 * hfs_constants::b * xi_prime(j, jp);
    auto hg0 = w6j(half, half, 0, n, f, jp) * w6j(half, half, 0, n, f, j);
    auto hg1 = w6j(half, half, 1, n, f, jp) * w6j(half, half, 1, n, f, j);

    return coeff * (hg1 - hg0);
}

static inline dcomplex hfs_tensor_angular(double j, double m_j, double jp,
    double m_jp, double f, double m_f,
    double fp, double m_fp, double n,
    double m_n, double m_i, double m_s) {
    const spin i = half, s = half;
    dcomplex w0 = w3j(i, jp, fp, m_i, m_jp, -m_fp);
    dcomplex w1 = w3j(s, n, jp, m_s, m_n, -m_jp);
    dcomplex w2 = w3j(i, j, f, m_i, m_j, -m_f);
    dcomplex w3 = w3j(s, n, j, m_s, m_n, -m_j);
    return w0 * w1 * w2 * w3 * m_i * m_s * parity(-m_j - m_jp);
}

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
    if (m_f != m_fp || f != fp) return 0;
    if (n != s2.n) return 0;

    const spin i = half, s = half;

    dcomplex prf = hfs_constants::c * xi_prime(jp, j);

    dcomplex htensor = 0;

    for (double m_j = -j; m_j <= j; m_j++) {
        for (double m_jp = -j; m_jp <= jp; m_jp++) {
            for (double m_n = -j; m_n <= n; m_n++) {
                // manually unrolled innermost loops
                htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, f, m_fp, n, m_n,
                    -half, -half);
                htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, f, m_fp, n, m_n,
                    -half, +half);
                htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, f, m_fp, n, m_n,
                    +half, -half);
                htensor += hfs_tensor_angular(j, m_j, jp, m_jp, f, m_f, f, m_fp, n, m_n,
                    +half, +half);
            }
        }
    }

    return prf * htensor;
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

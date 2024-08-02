#include <pch.h>
#include <aef/quantum.h>

aef::j_basis_vec aef::quantum::expectation_values_jsq(HyperfineCalculator& calc, int32_t E_idx) {
    //
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    j_basis_vec out;
#ifdef _WIN32
    SecureZeroMemory((void*)&out, sizeof(j_basis_vec));
#else
    //explicit_bzero((void*)&out, sizeof(j_basis_vec));
    memset((void*)&out, 0, sizeof(j_basis_vec));
#endif
    double prob_tot = 0;
    for (int32_t kidx = 0; kidx < calc.nBasisElts; kidx++) {
        const double prob = std::norm(state_vec[kidx]);

        if (prob < std::numeric_limits<double>::epsilon()) {
            continue;
        }

        prob_tot += prob;
        j_basis_vec bs_ket = calc.basis[kidx];
        // note that angular momenta
        out.n += prob * bs_ket.n * (bs_ket.n + 1);
        out.j += prob * bs_ket.j * (bs_ket.j + 1);
        out.f += prob * bs_ket.f * (bs_ket.f + 1);
        out.m_f += prob * bs_ket.m_f;
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }

    out.n = invert_qsq(out.n / prob_tot);
    out.j = invert_qsq(out.j / prob_tot);
    out.f = invert_qsq(out.f / prob_tot);
    out.m_f /= prob_tot;

    return out;
}

double aef::quantum::expect_parity(HyperfineCalculator& calc, int32_t E_idx) {
    double ex_parity = 0.0;
    double prob_tot = 0.0;
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    for (int32_t kidx = 0; kidx < calc.nBasisElts; kidx++) {
        const double prob = std::norm(state_vec[kidx]);

        if (prob < std::numeric_limits<double>::epsilon()) {
            continue;
        }

        prob_tot += prob;
        j_basis_vec bs_ket = calc.basis[kidx];
        ex_parity += prob * std::pow(-1, bs_ket.n);
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }
    return ex_parity / prob_tot;
}

#include <pch.h>
#include <aef/quantum.h>
#include <aef/MolecularSystem.h>

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

aef::universal_diatomic_basis_vec aef::quantum::expectation_values_jsq(aef::MolecularSystem& calc, int32_t E_idx) {
    //
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    auto* mcalc = calc.get_calc();
    aef::universal_diatomic_basis_vec out;
#ifdef _WIN32
    SecureZeroMemory((void*)&out, sizeof(j_basis_vec));
#else
    //explicit_bzero((void*)&out, sizeof(j_basis_vec));
    memset((void*)&out, 0, sizeof(j_basis_vec));
#endif
    out.type = mcalc->get_basis_ket(0).type;
    
    double prob_tot = 0;
    for (int32_t kidx = 0; kidx < calc.nBasisElts; kidx++) {
        const double prob = std::norm(state_vec[kidx]);

        if (prob < std::numeric_limits<double>::epsilon()) {
            continue;
        }

        prob_tot += prob;
        aef::universal_diatomic_basis_vec bs_ket = mcalc->get_basis_ket(kidx);
        // note that angular momenta
        out.i1 += prob * bs_ket.i1 * (bs_ket.i1 + 1);
        out.i2 += prob * bs_ket.i2 * (bs_ket.i2 + 1);
        out.s += prob * bs_ket.s * (bs_ket.s + 1);
        out.l += prob * bs_ket.l * (bs_ket.l + 1);
        out.r += prob * bs_ket.r * (bs_ket.r + 1);
        out.n += prob * bs_ket.n * (bs_ket.n + 1);
        out.j += prob * bs_ket.j * (bs_ket.j + 1);
        out.f_1 += prob * bs_ket.f_1 * (bs_ket.f_1 + 1);
        out.f += prob * bs_ket.f * (bs_ket.f + 1);
        out.m_f += prob * bs_ket.m_f;
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }

    // invert squaring
    out.i1 = invert_qsq(out.i1 / prob_tot);
    out.i2 = invert_qsq(out.i2 / prob_tot);
    out.l = invert_qsq(out.l / prob_tot);
    out.s = invert_qsq(out.s / prob_tot);
    out.r = invert_qsq(out.r / prob_tot);
    out.n = invert_qsq(out.n / prob_tot);
    out.j = invert_qsq(out.j / prob_tot);
    out.f_1 = invert_qsq(out.f_1 / prob_tot);
    out.f = invert_qsq(out.f / prob_tot);
    out.m_f /= prob_tot;

    return out;
}

double aef::quantum::expect_parity(aef::MolecularSystem& calc, int32_t E_idx) {
    double ex_parity = 0.0;
    double prob_tot = 0.0;
    // TODO come up with better solution
    auto* mcalc = calc.get_calc();
    Eigen::VectorXcd state_vec = calc.Vs.col(E_idx);
    for (int32_t kidx = 0; kidx < calc.nBasisElts; kidx++) {
        const double prob = std::norm(state_vec[kidx]);

        if (prob < std::numeric_limits<double>::epsilon()) {
            continue;
        }

        prob_tot += prob;
        aef::universal_diatomic_basis_vec bs_ket = mcalc->get_basis_ket(kidx);
        ex_parity += prob * std::pow(-1, bs_ket.n);
    }

    if (prob_tot > (1 + std::numeric_limits<double>::epsilon() * 100000)) {
        DebugBreak();
    }
    return ex_parity / prob_tot;
}

double aef::quantum::calculate_transition_rate(transition_type type, unsigned order, double energy, dcomplex mat_elt) {
    aef::quantum::transition_information tsn(type, order, energy, mat_elt);
    
    return tsn.A;
}

Eigen::MatrixXcd aef::orient_diag::makeOrientationDiagonalizer(aef::MolecularSystem& sys) {
    constexpr double inv_sqrt2 = std::numbers::sqrt2 / 2.0;
    Eigen::MatrixXcd orientEnergyMatrix;

    // convert to cartesian
    using namespace std::complex_literals;
    Eigen::MatrixXcd& dz = sys.d10;
    Eigen::MatrixXcd dx = (sys.d1t - sys.d11) * inv_sqrt2;
    Eigen::MatrixXcd dy = (sys.d1t + sys.d11) * 1i * inv_sqrt2;

    constexpr double E_dz = 40;
    constexpr double E_dx = 20;
    constexpr double E_dy = 10;

    orientEnergyMatrix = E_dz * dz + E_dx * dx + E_dy * dy;
    std::cout << fmt::format(
        "Orientation diagonalizer coeffs are E_dz = {} MHz, E_dx = {} MHz, E_dy = {} MHz",
        E_dz, E_dx, E_dy) << std::endl;
    return orientEnergyMatrix;
}
aef::ResultCode aef::orient_diag::diagonalize(aef::MolecularSystem& sys, Eigen::MatrixXcd& orientEnergyMatrix, Eigen::MatrixXcd* vals) {
    if (!vals) {
        vals = new Eigen::MatrixXcd();
        vals->resizeLike(sys.H_tot);
        vals->setZero();
    }


    sys.H_tot += sys.H_dev;
    *vals = sys.H_tot + orientEnergyMatrix;

    auto rc = aef::matrix::diagonalize(*vals, sys.Es, sys.Vs);
    assert("Diagonalization failed", aef::succeeded(rc));
    rc = aef::matrix::group_action(*vals, sys.Vs, sys.H_tot);
    assert("Eigenstate correction failed", aef::succeeded(rc));

    sys.Es = vals->diagonal();
    return rc;
}

aef::quantum::transition_information::transition_information(transition_type type_, unsigned order_, double freq_, dcomplex mat_elt_):
    type(type_), order(order_), mat_elt(mat_elt_), freq_MHz(freq_), calcs_done(false)
{
    A = B = t = f = std::nan("");
}

aef::ResultCode aef::quantum::transition_information::calculate() {
    A = base_rate() * std::norm(this->mat_elt);
    t = 1 / A;
    using namespace std::numbers;
    using namespace constants;
    constexpr auto debye = unit_conversion::C_m_per_D;
    double freq_Hz_cubed = freq_MHz * freq_MHz * freq_MHz * 1E18;
    constexpr auto F_v_coeff = 2 * h / (c*c);
    {
        constexpr auto A_1GHz = 1E-10;
        constexpr auto freq_1GHz = 1E9;
        constexpr auto freq_1GHz_cubed = freq_1GHz * freq_1GHz * freq_1GHz;
        constexpr auto F_v_1GHz = F_v_coeff * freq_1GHz_cubed;
        constexpr auto B_1GHz = A_1GHz / F_v_1GHz;
    }
    double F_v = 2 * h * freq_Hz_cubed;
    //B = A / F_v;
    {
        constexpr double B_coeff_si = 4 * pi * pi * pi / (3 * epsilon_naught * h * h);
        constexpr double B_coeff_D2 = B_coeff_si * debye * debye;
        B = B_coeff_D2 * std::norm(mat_elt);
    }


    return aef::ResultCode::Unimplemented;
}

double aef::quantum::transition_information::Energy_eV() const {
    double E_J = Energy_J();
    return E_J / constants::e;
}

double aef::quantum::transition_information::Energy_J() const {
    double f_Hz = freq_MHz * 1E6;
    return f_Hz * constants::h;
}

double aef::quantum::transition_information::wavelength_nm() const {
    double f_Hz = freq_MHz * 1E6;
    double l_m = constants::c / f_Hz;
    return l_m * 1E9;
}

double aef::quantum::transition_information::wavenumber_inv_cm() const {
    return freq_MHz /  unit_conversion::MHz_per_inv_cm;
}

double aef::quantum::transition_information::base_rate() const {
    assert("Higher order transitions not implemented yet", order == 1);

    using namespace std::numbers;
    using aef::quantum::transition_type::E;
    using aef::quantum::transition_type::M;
    using namespace constants;

    const double omega = 2 * pi * freq_MHz * 1E6; // rad/s
    const double k_rad = omega / c; // units: rad / m
    constexpr auto debye = unit_conversion::C_m_per_D;
    constexpr auto mub = constants::e * hbar / (2 * m_e);
    constexpr auto mmm = mub / h;

    if (order == 1 && type == E) {
        
        constexpr auto coeff = debye*debye *  pi / (3 * hbar * epsilon_naught);
        {
            constexpr double omega_700THz = 7E14 * 2 * pi;
            constexpr double k_rad_700THz = omega_700THz / c;
            constexpr double k_700THz_cubed = k_rad_700THz * k_rad_700THz * k_rad_700THz;
            constexpr double est_A_700THz = coeff * k_700THz_cubed;
        }

        {
            constexpr double omega_7GHz = 7E9 * 2 * pi;
            constexpr double k_rad_7GHz = omega_7GHz / c;
            constexpr double k_7GHz_cubed = k_rad_7GHz * k_rad_7GHz * k_rad_7GHz;
            constexpr double est_A_7GHz = coeff * k_7GHz_cubed;
        }

        return coeff *  k_rad * k_rad * k_rad;
    }

    if (order == 1 && type == M) {
        constexpr double J_per_T_from_MHz_per_T = h * 1E6;
        constexpr double coeff = J_per_T_from_MHz_per_T * J_per_T_from_MHz_per_T * mu_naught / (3 * pi * hbar);
        {
            constexpr double omega_700THz = 7E14 * 2 * pi;
            constexpr double k_rad_700THz = omega_700THz / c;
            constexpr double k_700THz_cubed = k_rad_700THz * k_rad_700THz * k_rad_700THz;
            constexpr double est_A_700THz = coeff * k_700THz_cubed *mu_bohr* mu_bohr;
        }
        return coeff * k_rad * k_rad * k_rad;
    }

    {
        constexpr double e1_m1_ratio = (debye* debye * pi / (3 * hbar * epsilon_naught)) / (mub * mub * mu_naught / (3 * pi * hbar));
        constexpr double norm_e1_m1_ratio = e1_m1_ratio * (alpha * alpha);
    }

    return 0.0;
}

double aef::quantum::transition_information::calc_A() const {
    const double base = base_rate();

    if (isnan(base)) {
        MessageBoxA(NULL, "FUK", "FUK base null", 0);
    }

    return base * std::norm(mat_elt);
}

double aef::quantum::transition_information::calc_f() const {
    assert("Higher order transitions not implemented yet", order == 1);

    using namespace std::numbers;
    using aef::quantum::transition_type::E;
    using aef::quantum::transition_type::M;
    using namespace constants;

    const double omega = 2 * pi * freq_MHz * 1E6; // rad/s
    const double k_rad = omega / c; // units: rad / m
    constexpr auto debye = unit_conversion::C_m_per_D;
    constexpr auto mub = constants::e * hbar / (2 * m_e);

    if (order == 1 && type == E) {
        // mat elt units: D, for JTS needs to be m
        constexpr double m_from_D_per_e = unit_conversion::C_m_per_D / constants::e;
        const double S = std::norm(m_from_D_per_e * mat_elt); // units: m^2

        constexpr double coeff = 2 * m_e / hbar; // units: kg / (J*s)
        {
            constexpr double omega_1GHz = 2 * pi * 1E9;
            constexpr double melt_1GHz = 1; // Debye
            constexpr double S_1GHz = std::norm(a0 * melt_1GHz);// / constants::e);
            constexpr double f_1GHz = coeff * S_1GHz * omega_1GHz;
        }

        return coeff * omega * S;
    }

    if (order == 1 && type == M) {

    }

    return 0.0;
}

auto fmt::formatter<aef::quantum::transition_information>::format(tsn_ifo tsn, format_context& ctx) const {
    double A = tsn.calc_A();
    const char* tsn_type = (tsn.type == aef::quantum::transition_type::E) ? "E" : "M";
    return fmt::formatter<std::string>::format(
        fmt::format("{}{} transition dE={}, {}", tsn_type, tsn.order, tsn.Energy_eV(), A), 
        ctx);
}
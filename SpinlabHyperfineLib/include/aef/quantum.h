#pragma once
#ifndef _AEF_QUANTUM_H
#define _AEF_QUANTUM_H 1

#include <aef/aef.h>
#include <aef/MolecularSystem.h>

namespace aef::quantum {
    /// <summary>
    /// Performs the quantum q*(q+1) squaring for angular momenta
    /// </summary>
    /// <param name="q"> q </param>
    /// <returns>q*(q+1)</returns>
    template <class T> T qsq(T q) {
        return q * (q + (T)1);
    }


    /// <summary>
    /// Inverts the quantum q*(q+1) squaring for angular momenta
    /// </summary>
    /// <param name="expect_qsq">the expectation value of &lt;q(q+1)&gt; </param>
    /// <returns>the effective expectation of q</returns>
    template <class T> T invert_qsq(T expect_qsq) {
        return (std::sqrt(4 * expect_qsq + 1.0) - 1.0) / 2.0;
    }



    /// <summary>
/// Calculates the expectation values of the basis operators.  The squared
/// </summary>
/// <param name="calc">HyperfineCalculator: contains operator matrix elements
/// and states</param> <param name="E_idx">the index of Energy level to
/// calculate with</param> <returns></returns>
    aef::j_basis_vec expectation_values_jsq(HyperfineCalculator& calc, int32_t E_idx);
    double expect_parity(HyperfineCalculator& calc, int32_t E_idx);

    aef::universal_diatomic_basis_vec expectation_values_jsq(aef::MolecularSystem &sys, int32_t E_idx);
    double expect_parity(aef::MolecularSystem& calc, int32_t E_idx);

    enum class transition_type {
        ELECTRIC,
        MAGNETIC,
        E = ELECTRIC,
        M = MAGNETIC
    };

    struct transition_information {
        transition_type type;
        unsigned order;

        double freq_MHz; // MHz, technically frequency f, not omega
        dcomplex mat_elt;
        bool calcs_done;
        
        // transition rates
        double A; // Einstein A coeff / decay rate, Hz
        double B; // Einstein B coeff, ??
        double f; // oscillator strength, dimensionless
        double t; // lifetime, = 1 / A

        /// <summary>
        /// 
        /// </summary>
        /// <param name="type_"></param>
        /// <param name="order_"></param>
        /// <param name="f_"></param>
        /// <param name="mat_elt_">The line matrix element with the following units:
        /// * for E1 transitions: Debye
        /// * for M1 transitions: MHz/T
        /// * other tranisitions not yet supported
        /// </param>
        transition_information(transition_type type_, unsigned order_, double f_, dcomplex mat_elt_);

        /// <summary>
        /// Calculate the transition rates, filling out A, B, f, and t
        /// </summary>
        /// <returns></returns>
        aef::ResultCode calculate();

        double Energy_eV() const;
        double Energy_J() const;
        double wavelength_nm() const;
        double wavenumber_inv_cm() const;
        double base_rate() const;
        double calc_A() const;
        double calc_f() const;
    };

    /// <summary>
    /// Calculates the transition rates
    /// </summary>
    /// <param name="type"></param>
    /// <param name="order"></param>
    /// <param name="energy"></param>
    /// <param name="mat_elt"></param>
    /// <returns></returns>
    double calculate_transition_rate(transition_type type, unsigned order, double energy, dcomplex mat_elt);
};

template <> struct fmt::formatter<aef::quantum::transition_information> : fmt::formatter<std::string> {
    using tsn_ifo = aef::quantum::transition_information;
    auto format(tsn_ifo tsn, format_context& ctx) const;
};

namespace aef::orient_diag {
    constexpr double E_dz = 40;
    constexpr double E_dx = 20;
    constexpr double E_dy = 10;


    Eigen::MatrixXcd makeOrientationDiagonalizer(aef::MolecularSystem& sys);
    aef::ResultCode diagonalize(aef::MolecularSystem& sys, Eigen::MatrixXcd &orientEnergyMatrix, Eigen::MatrixXcd *work=nullptr);
};

#endif //_AEF_QUANTUM_H
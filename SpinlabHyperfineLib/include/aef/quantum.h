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
};


#endif //_AEF_QUANTUM_H
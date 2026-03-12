/*
 * Copyright (c) 2008-2026  Nakata, Maho  All rights reserved.
 * BSD 2-Clause License (see mplapack_arithmetic_params.h for full text).
 */

#pragma once
#ifndef MPLAPACK_ARITHMETIC_PARAMS_GMP_H
#define MPLAPACK_ARITHMETIC_PARAMS_GMP_H

// Precondition: mplapack_arithmetic_params.h has been included.
// Precondition: <gmpxx.h> has been included via <mpblas.h>.
//   mpf_class, mpf_get_prec, mpf_mul_2exp, mpf_div_2exp must be in scope.
//
// GMP mpf does not have a fixed exponent range.  MPLAPACK uses a conservative
// safe exponent bound = max(mp_exp_t) / 2 to provide a stable range for
// the LAPACK underflow/overflow conventions.  All values are RUNTIME-DERIVED.

#include <limits>

namespace mplapack {
namespace detail {

    // Specialization of int_pow_base2 for mpf_class.
    // Uses mpf_mul_2exp / mpf_div_2exp for exact O(1) shifting.
    template <> inline mpf_class int_pow_base2<mpf_class>(arithmetic_int n) {
        const mpf_class one(1.0);
        mpf_class result;
        result.set_prec(mpf_get_prec(one.get_mpf_t()));
        if (n >= 0) {
            mpf_mul_2exp(result.get_mpf_t(), one.get_mpf_t(), static_cast<mp_bitcnt_t>(n));
        } else {
            mpf_div_2exp(result.get_mpf_t(), one.get_mpf_t(), static_cast<mp_bitcnt_t>(-n));
        }
        return result;
    }

} // namespace detail

namespace detail_gmp {

    // Conservative safe exponent bound for GMP mpf (identical to the historical
    // get_max_safe_exponent() in Rlamch.cpp).
    inline mp_bitcnt_t safe_emax() noexcept { return static_cast<mp_bitcnt_t>(std::numeric_limits<mp_exp_t>::max() / 2); }

} // namespace detail_gmp

// ---------------------------------------------------------------------------
// ArithmeticParams<mpf_class>
//
// GMP mpf arithmetic uses truncation (chopping), so:
//   eps  = 2^(-prec)    (rnd = 0.0)
//   emin = 1 - E        where E = safe_emax(), so rmin = 2^(-E)
//   emax = E
//   rmin = 2^(-E)
//   rmax = (1 - eps) * 2^E
//
// These semantics match the original RlamchXxx_gmp implementations.
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<mpf_class> get_arithmetic_params<mpf_class>() {
    ArithmeticParams<mpf_class> p;

    const mpf_class one(1.0);
    const mp_bitcnt_t nbits = mpf_get_prec(one.get_mpf_t());
    const mp_bitcnt_t E = detail_gmp::safe_emax();

    // eps = 2^(-nbits)
    p.eps.set_prec(nbits);
    mpf_div_2exp(p.eps.get_mpf_t(), one.get_mpf_t(), nbits);

    p.base = mpf_class(2.0);
    p.prec = p.eps * p.base;
    p.t.set_prec(nbits);
    mpf_set_ui(p.t.get_mpf_t(), nbits);
    p.rnd = mpf_class(0.0); // GMP truncates

    // emin = 1 - E, emax = E.
    // Use mpf_set_ui for exact integer conversion — casting through double would
    // round LONG_MAX/2 (which is 2^62-1, not a power of two) to 2^62.
    p.emax.set_prec(nbits);
    mpf_set_ui(p.emax.get_mpf_t(), E);
    p.emin.set_prec(nbits);
    mpf_sub(p.emin.get_mpf_t(), mpf_class(1.0).get_mpf_t(), p.emax.get_mpf_t());

    // rmin = 2^(-E)
    p.rmin.set_prec(nbits);
    mpf_div_2exp(p.rmin.get_mpf_t(), one.get_mpf_t(), E);

    // rmax = (1 - eps) * 2^E
    p.rmax.set_prec(nbits);
    const mpf_class mant = one - p.eps;
    mpf_mul_2exp(p.rmax.get_mpf_t(), mant.get_mpf_t(), E);

    // sfmin: netlib DLAMCH('S') logic
    p.sfmin = p.rmin;
    const mpf_class small = one / p.rmax;
    if (small >= p.sfmin)
        p.sfmin = small * (one + p.eps);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<mpf_class>
//
// All exponents computed at runtime from current precision and safe_emax().
// int_pow_base2<mpf_class> uses mpf_mul_2exp / mpf_div_2exp (exact, no
// std::pow with floating exponent).
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<mpf_class> get_blue_scaling_params<mpf_class>() {
    BlueScalingParams<mpf_class> q;

    const mpf_class one(1.0);
    const arithmetic_int nbits = static_cast<arithmetic_int>(mpf_get_prec(one.get_mpf_t()));
    const arithmetic_int E = static_cast<arithmetic_int>(detail_gmp::safe_emax());

    // GMP emin/emax in LAPACK convention:
    //   emin = 1 - E  =>  rmin = base^(emin-1) = 2^(-E)
    //   emax = E
    const arithmetic_int emin = arithmetic_int(1) - E;
    const arithmetic_int emax = E;

    q.exp_tsml = detail::ceildiv2(emin - 1);
    q.exp_tbig = detail::floordiv2(emax - nbits + 1);
    q.exp_ssml = -detail::floordiv2(emin - nbits);
    q.exp_sbig = -detail::ceildiv2(emax + nbits - 1);

    q.tsml = detail::int_pow_base2<mpf_class>(q.exp_tsml);
    q.tbig = detail::int_pow_base2<mpf_class>(q.exp_tbig);
    q.ssml = detail::int_pow_base2<mpf_class>(q.exp_ssml);
    q.sbig = detail::int_pow_base2<mpf_class>(q.exp_sbig);

    return q;
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_GMP_H

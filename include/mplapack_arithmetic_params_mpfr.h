/*
 * Copyright (c) 2008-2026  Nakata, Maho  All rights reserved.
 * BSD 2-Clause License (see mplapack_arithmetic_params.h for full text).
 */

#pragma once
#ifndef MPLAPACK_ARITHMETIC_PARAMS_MPFR_H
#define MPLAPACK_ARITHMETIC_PARAMS_MPFR_H

// Precondition: mplapack_arithmetic_params.h has been included.
// Precondition: mpreal.h (mpfrc++) has been included via <mpblas.h>.
//   The mpfr::mpreal type and mul_2si() must be in scope.
//
// All values are RUNTIME-DERIVED from the mpreal default precision and the
// current MPFR global exponent range.  This matches Rlamch_mpfr semantics.

namespace mplapack {
namespace detail {

    // Specialization of int_pow_base2 for mpfr::mpreal.
    // Uses mul_2si for exact, O(1) power-of-2 scaling at MPFR precision.
    template <> inline mpfr::mpreal int_pow_base2<mpfr::mpreal>(arithmetic_int n) {
        mpfr::mpreal one(1.0);
        return mul_2si(one, static_cast<mp_exp_t>(n));
    }

} // namespace detail

// ---------------------------------------------------------------------------
// ArithmeticParams<mpfr::mpreal>
//
// Every field is derived at call time from:
//   - mpreal::get_default_prec()  (current MPFR working precision)
//   - mpfr_get_emin() / mpfr_get_emax()
//   - mpfr_get_default_rounding_mode()
//
// IMPORTANT: The exponent range is interpreted in the MPLAPACK convention:
//   rmin = 2^(emin-1)  where emin = mpfr_get_emin().
//   rmax = (1 - eps) * 2^emax  where emax = mpfr_get_emax().
//
// The caller is responsible for applying the MPFR exponent clamping that
// MPLAPACK normally applies before calling Rlamch (e.g., truncating the
// exponent range to ±(prec * 64) to avoid bisection pathologies).
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<mpfr::mpreal> get_arithmetic_params<mpfr::mpreal>() {
    using REAL = mpfr::mpreal;
    ArithmeticParams<REAL> p;

    const REAL one(1.0);
    const mp_prec_t nbits = one.get_prec();
    const mp_exp_t emin = one.get_emin();
    const mp_exp_t emax = one.get_emax();

    // eps = 2^(-nbits)  (unit roundoff for round-to-nearest)
    p.eps = mul_2si(one, -static_cast<mp_exp_t>(nbits));
    p.base = REAL(2.0);
    p.prec = p.eps * p.base;
    p.t = REAL(static_cast<double>(nbits));
    p.rnd = (mpfr_get_default_rounding_mode() == MPFR_RNDZ) ? REAL(0.0) : REAL(1.0);
    p.emin = REAL(static_cast<double>(emin));
    p.emax = REAL(static_cast<double>(emax));

    // rmin = 2^(emin-1)
    p.rmin = mul_2si(one, emin - mp_exp_t(1));

    // rmax = (1 - eps) * 2^emax
    p.rmax = mul_2si(one - p.eps, emax);

    // sfmin: netlib DLAMCH('S') logic
    p.sfmin = p.rmin;
    const REAL small = one / p.rmax;
    if (small >= p.sfmin)
        p.sfmin = small * (one + p.eps);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<mpfr::mpreal>
//
// Exponents derived at runtime from the current MPFR precision/range.
// Formulas:
//   exp_tsml = ceil((emin - 1) / 2)
//   exp_tbig = floor((emax - digits + 1) / 2)
//   exp_ssml = -floor((emin - digits) / 2)
//   exp_sbig = -ceil((emax + digits - 1) / 2)
//
// int_pow_base2<mpreal> uses mul_2si (exact, no floating exponent).
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<mpfr::mpreal> get_blue_scaling_params<mpfr::mpreal>() {
    using REAL = mpfr::mpreal;
    BlueScalingParams<REAL> q;

    const REAL one(1.0);
    const arithmetic_int nbits = static_cast<arithmetic_int>(one.get_prec());
    const arithmetic_int emin = static_cast<arithmetic_int>(one.get_emin());
    const arithmetic_int emax = static_cast<arithmetic_int>(one.get_emax());

    q.exp_tsml = detail::ceildiv2(emin - 1);
    q.exp_tbig = detail::floordiv2(emax - nbits + 1);
    q.exp_ssml = -detail::floordiv2(emin - nbits);
    q.exp_sbig = -detail::ceildiv2(emax + nbits - 1);

    q.tsml = detail::int_pow_base2<REAL>(q.exp_tsml);
    q.tbig = detail::int_pow_base2<REAL>(q.exp_tbig);
    q.ssml = detail::int_pow_base2<REAL>(q.exp_ssml);
    q.sbig = detail::int_pow_base2<REAL>(q.exp_sbig);

    return q;
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_MPFR_H

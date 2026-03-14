/*
 * Copyright (c) 2008-2026  Nakata, Maho  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#pragma once
#ifndef MPLAPACK_ARITHMETIC_PARAMS_QD_H
#define MPLAPACK_ARITHMETIC_PARAMS_QD_H

// Precondition: mplapack_arithmetic_params.h must be included before this file.
// Precondition: qd_real type and its class constants (_eps, _min_normalized,
//   _max) must be in scope via <qd/qd_real.h> (pulled in through <mpblas.h>).
//
// qd_real is a quad-double type (unevaluated sum of 4 IEEE binary64 values).
// Its exponent range is bounded by the double exponent range. All values are
// derived from the QD library constants without any call to std::pow.
//
// Key parameters:
//   digits = std::numeric_limits<qd_real>::digits  (= 209)
//   emin   = frexp(qd_real::_min_normalized) exponent  (= -862, runtime)
//            NOTE: differs from std::numeric_limits<double>::min_exponent (= -1021).
//            qd_real::_min_normalized requires all four components to remain
//            within the normalized double range, which raises the effective
//            minimum exponent of the leading component.
//   emax   = std::numeric_limits<double>::max_exponent (= 1024, compile-time)
//
// Blue scaling exponents:
//   exp_tsml = ceil((-862-1)/2)         = -431
//   exp_tbig = floor((1024-209+1)/2)    =  408
//   exp_ssml = -floor((-862-209)/2)     =  536
//   exp_sbig = -ceil((1024+209-1)/2)    = -616

#include <cmath>
#include <limits>

namespace mplapack {
namespace detail {

    // Specialization of int_pow_base2 for qd_real.
    // The QD library provides ldexp(qd_real, int), which shifts a qd_real by an
    // exact power of 2.  All Blue exponents for qd_real are within [-1024, 1024],
    // well inside the double exponent range; no overflow/underflow can occur here.
    template <> inline qd_real int_pow_base2<qd_real>(arithmetic_int n) { return ::ldexp(qd_real(1.0), static_cast<int>(n)); }

} // namespace detail

// ---------------------------------------------------------------------------
// ArithmeticParams<qd_real>
// NOTE:
// The QD library upstream changed eps to "slightly more conservative values
// for eps." in NEWS (Changes for 2.2.4):
// https://sources.debian.org/src/qd/2.3.7-2.1/NEWS/
//
// For MPLAPACK Rlamch/compare stability, do not use qd_real::_eps here.
// Use fixed canonical E/P literals instead.
//
// sfmin still follows DLAMCH('S'):
//   sfmin = max(rmin, (1/rmax)*(1+eps))
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<qd_real> get_arithmetic_params<qd_real>() {
    ArithmeticParams<qd_real> p;

    const qd_real one(1.0);
    const qd_real two(2.0);

    // Canonical QD epsilon
    p.eps = qd_real(+0x1.fffffffffffffp-209,
                    +0x0.0000000000000p+0000,
                    +0x0.0000000000000p+0000,
                    +0x0.0000000000000p+0000);

    p.base = two;
    p.prec = p.eps * p.base;

    p.t = static_cast<arithmetic_int>(std::numeric_limits<qd_real>::digits);
    // QD uses IEEE double arithmetic internally; rounding occurs.
    p.rnd = one;

    int emin_int = 0;
    (void)std::frexp(qd_real::_min_normalized, &emin_int);
    // Do NOT assert against std::numeric_limits<qd_real>::min_exponent:
    // the QD library does not fully specialize std::numeric_limits<qd_real>,
    // and min_exponent may return 0 (the unspecialized default).
    // Use the frexp-derived value, consistent with the effective exponent
    // range of qd_real::_min_normalized.
    p.emin = static_cast<arithmetic_int>(emin_int);
    p.emax = static_cast<arithmetic_int>(std::numeric_limits<double>::max_exponent);

    p.rmin = qd_real::_min_normalized;
    p.rmax = qd_real::_max;

    // Rlamch("S"): safe minimum following netlib DLAMCH('S')
    p.sfmin = p.rmin;
    const qd_real small = one / p.rmax;
    if (small >= p.sfmin)
        p.sfmin = small * (one + p.eps);

    p.safmin = detail::compute_safmin<qd_real>(p.emin, p.emax);
    p.safmax = detail::compute_safmax<qd_real>(p.emin, p.emax);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<qd_real>
//
// All exponents are derived at runtime via the same ceildiv2/floordiv2
// helpers used for all other types.  ldexp is used for exact scaling.
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<qd_real> get_blue_scaling_params<qd_real>() {
    return make_blue_scaling_params(get_arithmetic_params<qd_real>());
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_QD_H

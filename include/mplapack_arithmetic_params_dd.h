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
#ifndef MPLAPACK_ARITHMETIC_PARAMS_DD_H
#define MPLAPACK_ARITHMETIC_PARAMS_DD_H

// Precondition: mplapack_arithmetic_params.h must be included before this file.
// Precondition: dd_real type and its class constants (_eps, _min_normalized,
//   _max) must be in scope via <qd/dd_real.h> (pulled in through <mpblas.h>).
//
// dd_real is a double-double type (unevaluated sum of 2 IEEE binary64 values).
// Its exponent range is bounded by the double exponent range. All values are
// derived from the QD library constants without any call to std::pow.
//
// Key parameters:
//   digits = std::numeric_limits<dd_real>::digits  (= 104)
//   emin   = frexp(dd_real::_min_normalized) exponent  (= -968, runtime)
//            NOTE: differs from std::numeric_limits<double>::min_exponent (= -1021).
//            dd_real::_min_normalized requires both components to remain
//            within the normalized double range, which raises the effective
//            minimum exponent of the leading component.
//   emax   = std::numeric_limits<double>::max_exponent (= 1024, compile-time)
//
// Blue scaling exponents:
//   exp_tsml = ceil((-968-1)/2)         = -484
//   exp_tbig = floor((1024-104+1)/2)    =  460
//   exp_ssml = -floor((-968-104)/2)     =  536
//   exp_sbig = -ceil((1024+104-1)/2)    = -564

#include <cmath>
#include <limits>

namespace mplapack {
namespace detail {

    // Specialization of int_pow_base2 for dd_real.
    // The QD library provides ldexp(dd_real, int) for exact binary scaling.
    // All Blue exponents for dd_real are within [-1024, 1024]; no
    // overflow/underflow can occur.
    template <> inline dd_real int_pow_base2<dd_real>(arithmetic_int n) { return ::ldexp(dd_real(1.0), static_cast<int>(n)); }

} // namespace detail

// ---------------------------------------------------------------------------
// ArithmeticParams<dd_real>
// NOTE:
// The QD library upstream changed eps to "slightly more conservative values
// for eps." in NEWS (Changes for 2.2.4):
// https://sources.debian.org/src/qd/2.3.7-2.1/NEWS/
//
// For MPLAPACK Rlamch/compare stability, do not use dd_real::_eps here.
// Use fixed canonical E/P literals instead.
//
// sfmin still follows DLAMCH('S'):
//   sfmin = max(rmin, (1/rmax)*(1+eps))
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<dd_real> get_arithmetic_params<dd_real>() {
    ArithmeticParams<dd_real> p;

    const dd_real one(1.0);
    const dd_real two(2.0);

    // Canonical DD epsilon for MPLAPACK
    p.eps = dd_real(+0x1.fffffffffffffp-105, +0x0.0000000000000p+0000);
    p.base = two;

    // PREC = EPS * BASE
    p.prec = p.eps * p.base;

    p.t = static_cast<arithmetic_int>(std::numeric_limits<dd_real>::digits);
    // DD uses IEEE double arithmetic internally; rounding occurs.
    p.rnd = one;

    int emin_int = 0;
    (void)std::frexp(dd_real::_min_normalized, &emin_int);
    // Do NOT assert against std::numeric_limits<dd_real>::min_exponent:
    // the QD library does not fully specialize std::numeric_limits<dd_real>
    // and min_exponent may return 0 (the unspecialized default).
    // Use the frexp-derived value, consistent with the effective exponent
    // range of dd_real::_min_normalized.
    p.emin = static_cast<arithmetic_int>(emin_int);
    p.emax = static_cast<arithmetic_int>(std::numeric_limits<double>::max_exponent);

    p.rmin = dd_real::_min_normalized;
    p.rmax = dd_real::_max;

    // Rlamch("S"): safe minimum following netlib DLAMCH('S')
    p.sfmin = p.rmin;
    const dd_real small = one / p.rmax;
    if (small >= p.sfmin)
        p.sfmin = small * (one + p.eps);

    p.safmin = detail::compute_safmin<dd_real>(p.emin, p.emax);
    p.safmax = detail::compute_safmax<dd_real>(p.emin, p.emax);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<dd_real>
//
// All exponents derived at runtime from the same sources used in
// get_arithmetic_params<dd_real>().
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<dd_real> get_blue_scaling_params<dd_real>() {
    return make_blue_scaling_params(get_arithmetic_params<dd_real>());
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_DD_H

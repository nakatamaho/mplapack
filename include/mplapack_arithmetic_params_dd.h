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
// Its exponent range is the same as IEEE binary64.  All values are derived
// from the QD library dd_real constants without any call to std::pow.
//
// Key parameters:
//   digits = std::numeric_limits<dd_real>::digits  (= 104)
//   emin   = frexp(dd_real::_min_normalized) exponent  (= -1021, runtime)
//   emax   = std::numeric_limits<double>::max_exponent (= 1024, compile-time)
//
// Blue scaling exponents (using ceildiv2/floordiv2 on arithmetic_int):
//   exp_tsml = ceil((-1021-1)/2)        = -511
//   exp_tbig = floor((1024-104+1)/2)    =  460
//   exp_ssml = -floor((-1021-104)/2)    =  563
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
//
// sfmin follows the same DLAMCH('S') logic as all other types:
//   sfmin = max(rmin, (1/rmax)*(1+eps))
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<dd_real> get_arithmetic_params<dd_real>() {
    ArithmeticParams<dd_real> p;

    const dd_real one(1.0);
    const dd_real two(2.0);

    p.eps = dd_real::_eps;
    p.base = two;
    p.prec = dd_real::_eps * two;
    p.t = dd_real(static_cast<double>(std::numeric_limits<dd_real>::digits));
    p.rnd = one; // DD uses IEEE double arithmetic internally; rounding occurs.

    // emin: use frexp of _min_normalized, matching legacy RlamchM_dd.
    int emin_int = 0;
    (void)std::frexp(dd_real::_min_normalized, &emin_int);
    p.emin = dd_real(static_cast<double>(emin_int));

    p.emax = dd_real(static_cast<double>(std::numeric_limits<double>::max_exponent));

    p.rmin = dd_real::_min_normalized;
    p.rmax = dd_real::_max;

    // sfmin: netlib DLAMCH('S') rule.
    p.sfmin = p.rmin;
    const dd_real small = one / p.rmax;
    if (small >= p.sfmin)
        p.sfmin = small * (one + p.eps);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<dd_real>
//
// All exponents derived at runtime from the same sources used in
// get_arithmetic_params<dd_real>().
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<dd_real> get_blue_scaling_params<dd_real>() {
    BlueScalingParams<dd_real> q;

    int emin_int = 0;
    (void)std::frexp(dd_real::_min_normalized, &emin_int);

    const arithmetic_int emin = static_cast<arithmetic_int>(emin_int);
    const arithmetic_int emax = static_cast<arithmetic_int>(std::numeric_limits<double>::max_exponent);
    const arithmetic_int digits = static_cast<arithmetic_int>(std::numeric_limits<dd_real>::digits);

    q.exp_tsml = detail::ceildiv2(emin - 1);
    q.exp_tbig = detail::floordiv2(emax - digits + 1);
    q.exp_ssml = -detail::floordiv2(emin - digits);
    q.exp_sbig = -detail::ceildiv2(emax + digits - 1);

    q.tsml = detail::int_pow_base2<dd_real>(q.exp_tsml);
    q.tbig = detail::int_pow_base2<dd_real>(q.exp_tbig);
    q.ssml = detail::int_pow_base2<dd_real>(q.exp_ssml);
    q.sbig = detail::int_pow_base2<dd_real>(q.exp_sbig);

    return q;
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_DD_H

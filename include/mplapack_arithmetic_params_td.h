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
#ifndef MPLAPACK_ARITHMETIC_PARAMS_TD_H
#define MPLAPACK_ARITHMETIC_PARAMS_TD_H

// Precondition: mplapack_arithmetic_params.h must be included before this file.
// Precondition: td_real type and its class constants (_eps, _min_normalized,
//   _max) must be in scope via <qd/td_real.h> (pulled in through <mpblas.h>).
//
// td_real is a triple-double type (unevaluated sum of 3 IEEE binary64 values,
// libQD3 1.6.0 or later).  Its exponent range is bounded by the double
// exponent range.
//
// Key parameters (libQD3 1.6.0):
//   digits = std::numeric_limits<td_real>::digits  (= 157)
//   emin   = frexp(td_real::_min_normalized) exponent  (= -915, runtime)
//            td_real::_min_normalized = 2^-916 keeps all three components in
//            the normalized double range, which raises the effective minimum
//            exponent of the leading component (dd: -968, qd: -862).
//   emax   = std::numeric_limits<double>::max_exponent (= 1024, compile-time)
//
// Blue scaling exponents:
//   exp_tsml = ceil((-915-1)/2)         = -458
//   exp_tbig = floor((1024-157+1)/2)    =  434
//   exp_ssml = -floor((-915-157)/2)     =  536
//   exp_sbig = -ceil((1024+157-1)/2)    = -590

#include <cmath>
#include <limits>

namespace mplapack {
namespace detail {

    // Specialization of int_pow_base2 for td_real.
    // libQD3 provides ldexp(td_real, int), which shifts a td_real by an exact
    // power of 2.  All Blue exponents for td_real are within [-1024, 1024].
    template <> inline td_real int_pow_base2<td_real>(arithmetic_int n) { return ::ldexp(td_real(1.0), static_cast<int>(n)); }

} // namespace detail

// ---------------------------------------------------------------------------
// ArithmeticParams<td_real>
// NOTE:
// dd and qd use fixed canonical E/P literals inherited from the decimal eps
// constants of the original QD library.  td has no such history: E is
// 2^-157 = 2^-digits, which is libQD3's td_real::_eps, written as a literal so
// Rlamch does not depend on the library constant.
//
// sfmin still follows DLAMCH('S'):
//   sfmin = max(rmin, (1/rmax)*(1+eps))
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<td_real> get_arithmetic_params<td_real>() {
    ArithmeticParams<td_real> p;

    const td_real one(1.0);
    const td_real two(2.0);

    p.eps = td_real(+0x1.0000000000000p-157, +0x0.0000000000000p+0000, +0x0.0000000000000p+0000);

    p.base = two;
    p.prec = p.eps * p.base;

    p.t = static_cast<arithmetic_int>(std::numeric_limits<td_real>::digits);
    // td_real uses IEEE double arithmetic internally; rounding occurs.
    p.rnd = one;

    int emin_int = 0;
    (void)std::frexp(td_real::_min_normalized, &emin_int);
    p.emin = static_cast<arithmetic_int>(emin_int);
    p.emax = static_cast<arithmetic_int>(std::numeric_limits<double>::max_exponent);

    p.rmin = td_real::_min_normalized;
    p.rmax = td_real::_max;

    // Rlamch("S"): safe minimum following netlib DLAMCH('S')
    p.sfmin = p.rmin;
    const td_real small = one / p.rmax;
    if (small >= p.sfmin)
        p.sfmin = small * (one + p.eps);

    p.safmin = detail::compute_safmin<td_real>(p.emin, p.emax);
    p.safmax = detail::compute_safmax<td_real>(p.emin, p.emax);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<td_real>
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<td_real> get_blue_scaling_params<td_real>() {
    return make_blue_scaling_params(get_arithmetic_params<td_real>());
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_TD_H

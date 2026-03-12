/*
 * Copyright (c) 2008-2026
 *      Nakata, Maho
 *      All rights reserved.
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
#ifndef MPLAPACK_ARITHMETIC_PARAMS_H
#define MPLAPACK_ARITHMETIC_PARAMS_H

// =============================================================================
// mplapack_arithmetic_params.h
//
// Common internal layer for mpblas and mplapack arithmetic parameter queries.
//
// Design goals:
//   - Single source of truth for floating-point format parameters.
//   - No dependency on Rlamch character-dispatch strings.
//   - No "machine" in naming; type-centric nomenclature throughout.
//   - MPFR and GMP values are always runtime-derived (no compile-time constants
//     for variable-precision types).
//   - Blue scaling exponents computed via integer arithmetic only (no
//     std::pow with floating exponent).
//   - Rlamch becomes a thin adapter over get_arithmetic_params<REAL>().
//   - No mpblas -> mplapack dependency introduced.
// =============================================================================

#include <cstdint>

namespace mplapack {

// ---------------------------------------------------------------------------
// arithmetic_int: integer type for exponent metadata.
// Deliberately NOT mplapackint to keep this layer independent of LAPACK
// integer width policies.
// ---------------------------------------------------------------------------
using arithmetic_int = std::intmax_t;

// ---------------------------------------------------------------------------
// Integer ceil/floor division by 2 (no floating-point involved).
// Used to compute Blue scaling exponents from emin/emax/digits.
// ---------------------------------------------------------------------------
namespace detail {

    constexpr arithmetic_int ceildiv2(arithmetic_int a) noexcept {
        // C++ truncates toward zero.
        // For a >= 0: ceil(a/2) = (a+1)/2.  For a < 0: ceil(a/2) = a/2.
        return (a >= 0) ? (a + 1) / 2 : a / 2;
    }

    constexpr arithmetic_int floordiv2(arithmetic_int a) noexcept {
        // For a >= 0: floor(a/2) = a/2.  For a < 0: floor(a/2) = (a-1)/2.
        return (a >= 0) ? a / 2 : (a - 1) / 2;
    }

    // Generic base-2 integer power via binary exponentiation.
    // Precondition: base is the type's representation of 2.
    // Never calls std::pow; uses only integer-counted multiplications.
    template <class REAL> REAL int_pow_base2(arithmetic_int n) {
        REAL result(1);
        if (n == 0)
            return result;
        if (n > 0) {
            REAL b(2);
            arithmetic_int e = n;
            while (e > 0) {
                if (e & arithmetic_int(1))
                    result = result * b;
                b = b * b;
                e >>= 1;
            }
            return result;
        }
        // n < 0: compute 2^(-n) then invert.
        REAL b(2);
        arithmetic_int e = -n;
        while (e > 0) {
            if (e & arithmetic_int(1))
                result = result * b;
            b = b * b;
            e >>= 1;
        }
        return REAL(1) / result;
    }

} // namespace detail

// ---------------------------------------------------------------------------
// ArithmeticParams<REAL>
//
// Mirrors the complete set of Rlamch selectors.
// Field correspondence (Rlamch selector -> field):
//   "E" -> eps    (unit roundoff, relative machine precision)
//   "S" -> sfmin  (safe minimum: 1/sfmin does not overflow)
//   "B" -> base   (radix; always 2 in MPLAPACK)
//   "P" -> prec   (precision = eps * base)
//   "N" -> t      (number of base-B mantissa digits)
//   "R" -> rnd    (rounding mode: 1.0 = rounds, 0.0 = chops)
//   "M" -> emin   (minimum exponent; rmin = base^(emin-1))
//   "U" -> rmin   (underflow threshold = minimum positive normal)
//   "L" -> emax   (largest exponent; rmax ~ base^emax)
//   "O" -> rmax   (overflow threshold = maximum finite value)
// ---------------------------------------------------------------------------
template <class REAL> struct ArithmeticParams {
    REAL eps;   // unit roundoff                [Rlamch "E"]
    REAL sfmin; // safe minimum                 [Rlamch "S"]
    REAL base;  // radix (= 2)                  [Rlamch "B"]
    REAL prec;  // eps * base                   [Rlamch "P"]
    REAL t;     // mantissa digit count         [Rlamch "N"]
    REAL rnd;   // rounding mode indicator      [Rlamch "R"]
    REAL emin;  // minimum exponent             [Rlamch "M"]
    REAL rmin;  // underflow threshold          [Rlamch "U"]
    REAL emax;  // largest exponent             [Rlamch "L"]
    REAL rmax;  // overflow threshold           [Rlamch "O"]
};

// ---------------------------------------------------------------------------
// BlueScalingParams<REAL>
//
// Scaling constants for the Blue algorithm (used in Rnrm2 / Clange etc.)
// All four exponents are stored as arithmetic_int for exact integer inspection.
//
// Formulas (base = 2, all exponents are integer-valued):
//   exp_tsml = ceil((emin - 1) / 2)
//   exp_tbig = floor((emax - digits + 1) / 2)
//   exp_ssml = -floor((emin - digits) / 2)
//   exp_sbig = -ceil((emax + digits - 1) / 2)
//
// The REAL values are base^exponent, computed via int_pow_base2 (no std::pow).
// ---------------------------------------------------------------------------
template <class REAL> struct BlueScalingParams {
    arithmetic_int exp_tsml; // exponent for tsml
    arithmetic_int exp_tbig; // exponent for tbig
    arithmetic_int exp_ssml; // exponent for ssml
    arithmetic_int exp_sbig; // exponent for sbig

    REAL tsml; // base^exp_tsml  (small-underflow threshold)
    REAL tbig; // base^exp_tbig  (large-overflow threshold)
    REAL ssml; // base^exp_ssml  (ssml^2 * x^2 safe near underflow)
    REAL sbig; // base^exp_sbig  (sbig^2 * x^2 safe near overflow)
};

// ---------------------------------------------------------------------------
// Primary template declarations.
// Explicit specializations live in the per-type headers below.
// ---------------------------------------------------------------------------
template <class REAL> ArithmeticParams<REAL> get_arithmetic_params();
template <class REAL> BlueScalingParams<REAL> get_blue_scaling_params();

} // namespace mplapack

// ---------------------------------------------------------------------------
// Pull in per-type specializations guarded by the active build target.
// Each header is self-contained and may only be included once per TU via
// the build-macro guard.
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_MPFR___)
#include "mplapack_arithmetic_params_mpfr.h"
#endif

#if defined(___MPLAPACK_BUILD_WITH_GMP___)
#include "mplapack_arithmetic_params_gmp.h"
#endif

#if defined(___MPLAPACK_BUILD_WITH_DOUBLE___)
#include "mplapack_arithmetic_params_double.h"
#endif

#if defined(___MPLAPACK_BUILD_WITH_BINARY80___)
#include "mplapack_arithmetic_params_binary80.h"
#endif

#if defined(___MPLAPACK_BUILD_WITH_BINARY128___)
#include "mplapack_arithmetic_params_binary128.h"
#endif

#if defined(___MPLAPACK_BUILD_WITH_QD___)
#include "mplapack_arithmetic_params_qd.h"
#endif

#if defined(___MPLAPACK_BUILD_WITH_DD___)
#include "mplapack_arithmetic_params_dd.h"
#endif

#endif // MPLAPACK_ARITHMETIC_PARAMS_H

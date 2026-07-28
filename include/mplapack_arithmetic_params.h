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
#include <cstdio>

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

    // Compute the internal scaling constant for the underflow side.
    // Matches the Fortran reference:
    //   safmin = radix ** max(minexponent-1, 1-maxexponent)
    // i.e. the larger (less negative) of the two exponent candidates is selected.
    // Note: this value is used internally for Blue scaling and is NOT
    // the underflow threshold returned by Rlamch("S").
    template <class REAL>
    inline REAL compute_safmin(arithmetic_int emin, arithmetic_int emax) {
        const arithmetic_int exp_from_emin = emin - arithmetic_int(1);
        const arithmetic_int exp_from_emax = arithmetic_int(1) - emax;

        const arithmetic_int exp_safmin =
            (exp_from_emin > exp_from_emax) ? exp_from_emin : exp_from_emax;

        return int_pow_base2<REAL>(exp_safmin);
    }

    // Compute the internal scaling constant for the overflow side.
    // Matches the Fortran reference:
    //   safmax = radix ** max(1-minexponent, maxexponent-1)
    // i.e. the larger of the two exponent candidates is selected.
    // Note: this value is used internally for Blue scaling and is NOT
    // the overflow threshold returned by Rlamch("O").
    template <class REAL>
    inline REAL compute_safmax(arithmetic_int emin, arithmetic_int emax) {
        const arithmetic_int exp_up_from_emin = arithmetic_int(1) - emin;
        const arithmetic_int exp_down_from_emax = emax - arithmetic_int(1);

        const arithmetic_int exp_safmax =
            (exp_up_from_emin > exp_down_from_emax) ? exp_up_from_emin : exp_down_from_emax;
        return int_pow_base2<REAL>(exp_safmax);
    }

    // GMP-specific exact power-of-two helper using mpf_mul_2exp/mpf_div_2exp.
#if defined(__GNU_MPXX_H__) || defined(__GMP_PLUSPLUS__) || defined(GMPXX_MKII_H)
    inline mpf_class int_pow_base2_gmp(arithmetic_int exp) {
        mpf_class x(1, gmpxx::default_mpf_precision_bits());
        if (exp >= 0) {
            mpf_mul_2exp(x.get_mpf_t(), x.get_mpf_t(), static_cast<mp_bitcnt_t>(exp));
        } else {
            mpf_div_2exp(x.get_mpf_t(), x.get_mpf_t(), static_cast<mp_bitcnt_t>(-exp));
        }
        return x;
    }

    // GMP-specific implementation path for the internal underflow-side scaling
    // constant.  Use the same exponent-side selection rule as the generic helper,
    // but construct the exact power-of-two value with mpf_mul_2exp/mpf_div_2exp.
    template <>
    inline mpf_class compute_safmin<mpf_class>(arithmetic_int emin, arithmetic_int emax) {
        const arithmetic_int exp_from_emin = emin - arithmetic_int(1);
        const arithmetic_int exp_from_emax = arithmetic_int(1) - emax;

        const arithmetic_int exp_safmin =
            (exp_from_emin < exp_from_emax) ? exp_from_emin : exp_from_emax;

        return int_pow_base2_gmp(exp_safmin);
    }

    // GMP special-case:
    // mpf_get_emin()/mpf_get_emax() may describe a much wider exponent interval
    // than the finite power-of-two scaling range we want for MPLAPACK's internal
    // safe-scaling constants.  For the overflow-side helper, using the larger
    // exponent distance can make the internal scaling constant unnecessarily
    // aggressive, so for mpf_class we conservatively use the smaller exponent-side
    // distance.
    template <>
    inline mpf_class compute_safmax<mpf_class>(arithmetic_int emin, arithmetic_int emax) {
        const arithmetic_int exp_up_from_emin = arithmetic_int(1) - emin;
        const arithmetic_int exp_down_from_emax = emax - arithmetic_int(1);

        const arithmetic_int exp_safmax =
            (exp_up_from_emin < exp_down_from_emax) ? exp_up_from_emin : exp_down_from_emax;

        return int_pow_base2_gmp(exp_safmax);
    }
#endif // defined(__GNU_MPXX_H__) || defined(__GMP_PLUSPLUS__) || defined(GMPXX_MKII_H)

    // Warn when exponent-range conditions required by the internal
    // Blue scaling guards are violated. For the derived internal
    // power-of-two guards, the intended ordering is safmin <= 1 <= safmax.
    // This check is mainly relevant to runtime-configurable exponent ranges
    // (for example, MPFR). Fixed-parameter backends normally keep emin/emax
    // constant and should not hit this path in normal use.
    inline void warn_if_invalid_blue_scaling_exp_range(arithmetic_int emin, arithmetic_int emax) {
        const bool invalid_range = (emin >= emax);
        const bool emin_above_one = (emin > 1);
        const bool emax_below_one = (emax < 1);

        if (!invalid_range && !emin_above_one && !emax_below_one)
            return;

        if (invalid_range) {
            std::fprintf(stderr,
			 "[mplapack WARNING] Blue scaling: invalid exponent range: emin=%jd, emax=%jd "
                         "(expected emin < emax); cannot derive consistent internal Blue-scaling thresholds.\n",
			 static_cast<std::intmax_t>(emin), static_cast<std::intmax_t>(emax));
            return;
        }

        if (emin_above_one) {
            std::fprintf(stderr,
			 "[mplapack WARNING] Blue scaling: emin=%jd, emax=%jd: emin > 1 violates the expected "
                         "internal Blue-scaling ordering safmin <= 1 <= safmax; the derived internal guards become "
                         "safmin > 1.0 and safmax < 1.0.\n",
                         static_cast<std::intmax_t>(emin), static_cast<std::intmax_t>(emax));

        }

        if (emax_below_one) {
            std::fprintf(stderr,
			 "[mplapack WARNING] Blue scaling: emin=%jd, emax=%jd: emax < 1 violates the expected "
                         "internal Blue-scaling ordering safmin <= 1 <= safmax; the derived internal guards become "
                         "safmin > 1.0 and safmax < 1.0.\n",
                         static_cast<std::intmax_t>(emin), static_cast<std::intmax_t>(emax));

        }
    }

    // Default conversion path for small integer Rlamch metadata.
    // This uses a double intermediate and is intended for backends where
    // t/emin/emax are well within the exactly representable range of double.
    template <class REAL> inline REAL to_rlamch_real_impl(arithmetic_int x, REAL *) { return REAL(static_cast<double>(x)); }

#if defined(MPLAPACK_BUILD_WITH_GMP)
    inline mpf_class to_rlamch_real_impl(arithmetic_int x, mpf_class *) {
        const std::string s = std::to_string(x);
        return mpf_class(s.c_str());
    }
#endif

#if defined(MPLAPACK_BUILD_WITH_MPFR)
    inline mpfrxx::mpfr_class to_rlamch_real_impl(arithmetic_int x, mpfrxx::mpfr_class *) {
        const std::string s = std::to_string(x);
        return mpfrxx::mpfr_class(s.c_str());
    }
#endif

    template <class REAL> inline REAL to_rlamch_real(arithmetic_int x) { return to_rlamch_real_impl(x, static_cast<REAL *>(nullptr)); }
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
    REAL rnd;   // rounding mode indicator      [Rlamch "R"]
    REAL rmin;  // underflow threshold          [Rlamch "U"]
    REAL rmax;  // overflow threshold           [Rlamch "O"]

    arithmetic_int t;    // mantissa digit count [Rlamch "N"]
    arithmetic_int emin; // minimum exponent     [Rlamch "M"]
    arithmetic_int emax; // largest exponent     [Rlamch "L"]

    REAL safmin; // internal underflow-side power-of-two scaling constant
    REAL safmax; // internal overflow-side power-of-two scaling constant
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

    REAL tsml; // threshold below which values are classified as small
    REAL tbig; // threshold above which values are classified as big
    REAL ssml; // scaling factor applied to small values
    REAL sbig; // scaling factor applied to big values
};

// Assumes ap.base == 2 and uses exact powers of two via int_pow_base2.
template <class REAL>
inline BlueScalingParams<REAL> make_blue_scaling_params(const ArithmeticParams<REAL> &ap) {
    detail::warn_if_invalid_blue_scaling_exp_range(ap.emin, ap.emax);
    BlueScalingParams<REAL> q;

    q.exp_tsml = detail::ceildiv2(ap.emin - arithmetic_int(1));
    q.exp_tbig = detail::floordiv2(ap.emax - ap.t + arithmetic_int(1));
    q.exp_ssml = -detail::floordiv2(ap.emin - ap.t);
    q.exp_sbig = -detail::ceildiv2(ap.emax + ap.t - arithmetic_int(1));

    q.tsml = detail::int_pow_base2<REAL>(q.exp_tsml);
    q.tbig = detail::int_pow_base2<REAL>(q.exp_tbig);
    q.ssml = detail::int_pow_base2<REAL>(q.exp_ssml);
    q.sbig = detail::int_pow_base2<REAL>(q.exp_sbig);
    return q;
}

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
#if defined(MPLAPACK_BUILD_WITH_MPFR)
#include "mplapack_arithmetic_params_mpfr.h"
#endif

#if defined(MPLAPACK_BUILD_WITH_GMP)
#include "mplapack_arithmetic_params_gmp.h"
#endif

#if defined(MPLAPACK_BUILD_WITH_DOUBLE)
#include "mplapack_arithmetic_params_double.h"
#endif

#if defined(MPLAPACK_BUILD_WITH_BINARY80)
#include "mplapack_arithmetic_params_binary80.h"
#endif

#if defined(MPLAPACK_BUILD_WITH_BINARY128)
#include "mplapack_arithmetic_params_binary128.h"
#endif

#if defined(MPLAPACK_BUILD_WITH_QD)
#include "mplapack_arithmetic_params_qd.h"
#endif

#if defined(MPLAPACK_BUILD_WITH_DD)
#include "mplapack_arithmetic_params_dd.h"
#endif

#endif // MPLAPACK_ARITHMETIC_PARAMS_H

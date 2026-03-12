/*
 * Copyright (c) 2008-2026  Nakata, Maho  All rights reserved.
 * BSD 2-Clause License (see mplapack_arithmetic_params.h for full text).
 */

#pragma once
#ifndef MPLAPACK_ARITHMETIC_PARAMS_BINARY80_H
#define MPLAPACK_ARITHMETIC_PARAMS_BINARY80_H

// Precondition: mplapack_arithmetic_params.h has been included.
// Precondition: mplapack_config.h (or equivalent) defining
//   MPLAPACK_BINARY80_MODE, mplapack_binary80_t, and the
//   MPLAPACK_BINARY80_MODE_* constants has been included via <mpblas.h>.

#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
#include <cfloat>

namespace mplapack {
namespace detail {

    // Specialization of int_pow_base2 for mplapack_binary80_t.
    // Uses scalbnl (standard C99/C++11) for exact power-of-2 scaling when
    // the exponent is within the hardware range; falls back to binary
    // exponentiation for very large exponents such as Blue scaling values.
    // For binary80 the exponent range [-16382, 16383] is so wide that the
    // Blue exponents always fit, so scalbnl is always safe here.
    template <> inline mplapack_binary80_t int_pow_base2<mplapack_binary80_t>(arithmetic_int n) {
        // scalbnl is exact for powers of 2 within the representable range.
        // arithmetic_int is std::intmax_t; safe to narrow since binary80 exponent
        // fits comfortably inside int.
        return scalbnl((mplapack_binary80_t)1.0L, static_cast<int>(n));
    }

} // namespace detail

// ---------------------------------------------------------------------------
// Macro helpers: resolve mant_dig / min_exp / max_exp / min / max for the
// active binary80 implementation mode.
// ---------------------------------------------------------------------------
namespace detail_b80 {

    // Returns mantissa digit count (= 64 for genuine binary80).
    inline constexpr int mant_dig() noexcept {
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
#if defined(FLT64X_MANT_DIG)
        return FLT64X_MANT_DIG;
#elif defined(__FLT64X_MANT_DIG__)
        return __FLT64X_MANT_DIG__;
#else
#error "FLT64X_MANT_DIG unavailable for MPLAPACK_BINARY80_MODE_FLOAT64X"
#endif
#elif MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
        return LDBL_MANT_DIG;
#else
#error "Unknown or disabled MPLAPACK_BINARY80_MODE"
#endif
    }

    inline constexpr int min_exp() noexcept {
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
#if defined(FLT64X_MIN_EXP)
        return FLT64X_MIN_EXP;
#elif defined(__FLT64X_MIN_EXP__)
        return __FLT64X_MIN_EXP__;
#else
#error "FLT64X_MIN_EXP unavailable"
#endif
#elif MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
        return LDBL_MIN_EXP;
#else
#error "Unknown or disabled MPLAPACK_BINARY80_MODE"
#endif
    }

    inline constexpr int max_exp() noexcept {
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
#if defined(FLT64X_MAX_EXP)
        return FLT64X_MAX_EXP;
#elif defined(__FLT64X_MAX_EXP__)
        return __FLT64X_MAX_EXP__;
#else
#error "FLT64X_MAX_EXP unavailable"
#endif
#elif MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
        return LDBL_MAX_EXP;
#else
#error "Unknown or disabled MPLAPACK_BINARY80_MODE"
#endif
    }

    inline mplapack_binary80_t rmin_val() noexcept {
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
#if defined(FLT64X_MIN)
        return (mplapack_binary80_t)FLT64X_MIN;
#elif defined(__FLT64X_MIN__)
        return (mplapack_binary80_t)__FLT64X_MIN__;
#else
#error "FLT64X_MIN unavailable"
#endif
#elif MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
        return (mplapack_binary80_t)LDBL_MIN;
#else
#error "Unknown or disabled MPLAPACK_BINARY80_MODE"
#endif
    }

    inline mplapack_binary80_t rmax_val() noexcept {
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
#if defined(FLT64X_MAX)
        return (mplapack_binary80_t)FLT64X_MAX;
#elif defined(__FLT64X_MAX__)
        return (mplapack_binary80_t)__FLT64X_MAX__;
#else
#error "FLT64X_MAX unavailable"
#endif
#elif MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
        return (mplapack_binary80_t)LDBL_MAX;
#else
#error "Unknown or disabled MPLAPACK_BINARY80_MODE"
#endif
    }

} // namespace detail_b80

// ---------------------------------------------------------------------------
// ArithmeticParams<mplapack_binary80_t>
//
// digits = 64, emin = -16381, emax = 16384  (x87 80-bit extended precision)
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<mplapack_binary80_t> get_arithmetic_params<mplapack_binary80_t>() {
    using T = mplapack_binary80_t;
    ArithmeticParams<T> p;

    const int kMantDig = detail_b80::mant_dig();

    // eps = 2^(-kMantDig)  (unit roundoff, round-to-nearest)
    T eps = (T)1.0L;
    for (int i = 0; i < kMantDig; ++i)
        eps /= (T)2.0L;

    const T rmin = detail_b80::rmin_val();
    const T rmax = detail_b80::rmax_val();

    p.eps = eps;
    p.base = (T)2.0L;
    p.prec = eps * (T)2.0L;
    p.t = (T)kMantDig;
    p.rnd = (T)1.0L; // IEEE rounds
    p.emin = (T)detail_b80::min_exp();
    p.rmin = rmin;
    p.emax = (T)detail_b80::max_exp();
    p.rmax = rmax;

    p.sfmin = rmin;
    const T small = (T)1.0L / rmax;
    if (small >= p.sfmin)
        p.sfmin = small * ((T)1.0L + eps);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<mplapack_binary80_t>
//
// emin = -16381, emax = 16384, digits = 64
//   exp_tsml = ceil((-16381-1)/2)      = -8191
//   exp_tbig = floor((16384-64+1)/2)   =  8160  (floor(16321/2))
//   exp_ssml = -floor((-16381-64)/2)   =  8223  (-floor(-8222.5) = 8223)
//   exp_sbig = -ceil((16384+64-1)/2)   = -8224  (-ceil(8223.5)   = -8224)
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<mplapack_binary80_t> get_blue_scaling_params<mplapack_binary80_t>() {
    BlueScalingParams<mplapack_binary80_t> q;

    const arithmetic_int emin = static_cast<arithmetic_int>(detail_b80::min_exp());
    const arithmetic_int emax = static_cast<arithmetic_int>(detail_b80::max_exp());
    const arithmetic_int digits = static_cast<arithmetic_int>(detail_b80::mant_dig());

    q.exp_tsml = detail::ceildiv2(emin - 1);
    q.exp_tbig = detail::floordiv2(emax - digits + 1);
    q.exp_ssml = -detail::floordiv2(emin - digits);
    q.exp_sbig = -detail::ceildiv2(emax + digits - 1);

    q.tsml = detail::int_pow_base2<mplapack_binary80_t>(q.exp_tsml);
    q.tbig = detail::int_pow_base2<mplapack_binary80_t>(q.exp_tbig);
    q.ssml = detail::int_pow_base2<mplapack_binary80_t>(q.exp_ssml);
    q.sbig = detail::int_pow_base2<mplapack_binary80_t>(q.exp_sbig);

    return q;
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_BINARY80_H

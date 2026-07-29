/*
 * Copyright (c) 2008-2026  Nakata, Maho  All rights reserved.
 * BSD 2-Clause License (see mplapack_arithmetic_params.h for full text).
 */

#pragma once
#ifndef MPLAPACK_ARITHMETIC_PARAMS_BINARY128_H
#define MPLAPACK_ARITHMETIC_PARAMS_BINARY128_H

// Precondition: mplapack_arithmetic_params.h has been included.
// Precondition: mplapack_config.h (or equivalent) defining
//   MPLAPACK_BINARY128_MODE, mplapack_binary128_t, and the
//   MPLAPACK_BINARY128_MODE_* constants has been included via <mpblas.h>.

#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#include <quadmath.h>
#endif
#include <cfloat>

namespace mplapack {
namespace detail_b128 {

    // ---------------------------------------------------------------------------
    // Mode-dispatched compile-time constants for binary128.
    // ---------------------------------------------------------------------------

    inline constexpr int mant_dig() noexcept {
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
#if defined(__FLT128_MANT_DIG__)
        return __FLT128_MANT_DIG__;
#elif defined(FLT128_MANT_DIG)
        return FLT128_MANT_DIG;
#else
#error "FLT128_MANT_DIG required for MPLAPACK_BINARY128_MODE_FLOAT128"
#endif
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
        return LDBL_MANT_DIG;
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#if defined(__FLT128_MANT_DIG__)
        return __FLT128_MANT_DIG__;
#elif defined(FLT128_MANT_DIG)
        return FLT128_MANT_DIG;
#else
#error "FLT128_MANT_DIG unavailable for MPLAPACK_BINARY128_MODE_QUADMATH"
#endif
#else
#error "Unknown or disabled MPLAPACK_BINARY128_MODE"
#endif
    }

    inline constexpr int min_exp() noexcept {
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
#if defined(__FLT128_MIN_EXP__)
        return __FLT128_MIN_EXP__;
#elif defined(FLT128_MIN_EXP)
        return FLT128_MIN_EXP;
#else
#error "FLT128_MIN_EXP required"
#endif
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
        return LDBL_MIN_EXP;
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#if defined(__FLT128_MIN_EXP__)
        return __FLT128_MIN_EXP__;
#elif defined(FLT128_MIN_EXP)
        return FLT128_MIN_EXP;
#else
#error "FLT128_MIN_EXP unavailable"
#endif
#else
#error "Unknown or disabled MPLAPACK_BINARY128_MODE"
#endif
    }

    inline constexpr int max_exp() noexcept {
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
#if defined(__FLT128_MAX_EXP__)
        return __FLT128_MAX_EXP__;
#elif defined(FLT128_MAX_EXP)
        return FLT128_MAX_EXP;
#else
#error "FLT128_MAX_EXP required"
#endif
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
        return LDBL_MAX_EXP;
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#if defined(__FLT128_MAX_EXP__)
        return __FLT128_MAX_EXP__;
#elif defined(FLT128_MAX_EXP)
        return FLT128_MAX_EXP;
#else
#error "FLT128_MAX_EXP unavailable"
#endif
#else
#error "Unknown or disabled MPLAPACK_BINARY128_MODE"
#endif
    }

    inline mplapack_binary128_t rmin_val() noexcept {
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
#if defined(__FLT128_MIN__)
        return (mplapack_binary128_t)__FLT128_MIN__;
#elif defined(FLT128_MIN)
        return (mplapack_binary128_t)FLT128_MIN;
#else
#error "FLT128_MIN required"
#endif
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
        return (mplapack_binary128_t)LDBL_MIN;
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#if defined(FLT128_MIN)
        return (mplapack_binary128_t)FLT128_MIN;
#elif defined(__FLT128_MIN__)
        return (mplapack_binary128_t)__FLT128_MIN__;
#else
#error "FLT128_MIN unavailable"
#endif
#else
#error "Unknown or disabled MPLAPACK_BINARY128_MODE"
#endif
    }

    inline mplapack_binary128_t rmax_val() noexcept {
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
#if defined(__FLT128_MAX__)
        return (mplapack_binary128_t)__FLT128_MAX__;
#elif defined(FLT128_MAX)
        return (mplapack_binary128_t)FLT128_MAX;
#else
#error "FLT128_MAX required"
#endif
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
        return (mplapack_binary128_t)LDBL_MAX;
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#if defined(FLT128_MAX)
        return (mplapack_binary128_t)FLT128_MAX;
#elif defined(__FLT128_MAX__)
        return (mplapack_binary128_t)__FLT128_MAX__;
#else
#error "FLT128_MAX unavailable"
#endif
#else
#error "Unknown or disabled MPLAPACK_BINARY128_MODE"
#endif
    }

} // namespace detail_b128

// ---------------------------------------------------------------------------
// ArithmeticParams<mplapack_binary128_t>
//
// digits = 113, emin = -16381, emax = 16384  (IEEE 754 binary128)
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<mplapack_binary128_t> get_arithmetic_params<mplapack_binary128_t>() {
    using T = mplapack_binary128_t;
    ArithmeticParams<T> p;

    const int kMantDig = detail_b128::mant_dig();

    // eps = 2^(-kMantDig)
    T eps = (T)1.0;
    for (int i = 0; i < kMantDig; ++i)
        eps /= (T)2.0;

    const T rmin = detail_b128::rmin_val();
    const T rmax = detail_b128::rmax_val();

    p.eps = eps;
    p.base = (T)2.0;
    p.prec = eps * (T)2.0;
    p.t = static_cast<arithmetic_int>(kMantDig);
    p.rnd = (T)1.0; // IEEE rounds
    p.emin = static_cast<arithmetic_int>(detail_b128::min_exp());
    p.rmin = rmin;
    p.emax = static_cast<arithmetic_int>(detail_b128::max_exp());
    p.rmax = rmax;

    p.sfmin = rmin;
    const T small = (T)1.0 / rmax;
    if (small >= p.sfmin)
        p.sfmin = small * ((T)1.0 + eps);

    p.safmin = detail::compute_safmin<T>(p.emin, p.emax);
    p.safmax = detail::compute_safmax<T>(p.emin, p.emax);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<mplapack_binary128_t>
//
// Reference exponents for IEEE binary128:
//   emin = -16381, emax = 16384, digits = 113
//   exp_tsml = ceil((-16381-1)/2)       = -8191
//   exp_tbig = floor((16384-113+1)/2)   =  8136
//   exp_ssml = -floor((-16381-113)/2)   =  8247
//   exp_sbig = -ceil((16384+113-1)/2)   = -8248
//
// Values are derived via make_blue_scaling_params(get_arithmetic_params<T>()).
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<mplapack_binary128_t> get_blue_scaling_params<mplapack_binary128_t>() {
    return make_blue_scaling_params(get_arithmetic_params<mplapack_binary128_t>());
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_BINARY128_H





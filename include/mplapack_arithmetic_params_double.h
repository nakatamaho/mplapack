/*
 * Copyright (c) 2008-2026  Nakata, Maho  All rights reserved.
 * BSD 2-Clause License (see mplapack_arithmetic_params.h for full text).
 */

#pragma once
#ifndef MPLAPACK_ARITHMETIC_PARAMS_DOUBLE_H
#define MPLAPACK_ARITHMETIC_PARAMS_DOUBLE_H

// Precondition: mplapack_arithmetic_params.h has been included.
// Precondition: <cfloat> is available (DBL_* macros).

#include <cfloat>

namespace mplapack {

// ---------------------------------------------------------------------------
// ArithmeticParams<double>
//
// All values derived from <cfloat> macros; no runtime computation required.
//   digits = DBL_MANT_DIG = 53
//   emin   = DBL_MIN_EXP  = -1021
//   emax   = DBL_MAX_EXP  = 1024
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<double> get_arithmetic_params<double>() {
    ArithmeticParams<double> p;

    // eps = 2^(-DBL_MANT_DIG)  (unit roundoff, round-to-nearest)
    double eps = 1.0;
    for (int i = 0; i < DBL_MANT_DIG; ++i)
        eps /= 2.0;

    p.eps = eps;
    p.base = 2.0;
    p.prec = eps * 2.0;
    p.t = static_cast<double>(DBL_MANT_DIG);
    p.rnd = 1.0; // IEEE 754 rounds
    p.emin = static_cast<double>(DBL_MIN_EXP);
    p.rmin = DBL_MIN;
    p.emax = static_cast<double>(DBL_MAX_EXP);
    p.rmax = DBL_MAX;

    // sfmin: netlib DLAMCH('S') logic
    p.sfmin = DBL_MIN;
    const double small = 1.0 / DBL_MAX;
    if (small >= p.sfmin)
        p.sfmin = small * (1.0 + eps);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<double>
//
// emin = -1021, emax = 1024, digits = 53
//   exp_tsml = ceil((-1021-1)/2)        = -511
//   exp_tbig = floor((1024-53+1)/2)     =  486
//   exp_ssml = -floor((-1021-53)/2)     =  537
//   exp_sbig = -ceil((1024+53-1)/2)     = -538
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<double> get_blue_scaling_params<double>() {
    BlueScalingParams<double> q;

    constexpr arithmetic_int emin = static_cast<arithmetic_int>(DBL_MIN_EXP);
    constexpr arithmetic_int emax = static_cast<arithmetic_int>(DBL_MAX_EXP);
    constexpr arithmetic_int digits = static_cast<arithmetic_int>(DBL_MANT_DIG);

    q.exp_tsml = detail::ceildiv2(emin - 1);
    q.exp_tbig = detail::floordiv2(emax - digits + 1);
    q.exp_ssml = -detail::floordiv2(emin - digits);
    q.exp_sbig = -detail::ceildiv2(emax + digits - 1);

    q.tsml = detail::int_pow_base2<double>(q.exp_tsml);
    q.tbig = detail::int_pow_base2<double>(q.exp_tbig);
    q.ssml = detail::int_pow_base2<double>(q.exp_ssml);
    q.sbig = detail::int_pow_base2<double>(q.exp_sbig);

    return q;
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_DOUBLE_H

/*
 * Copyright (c) 2008-2026  Nakata, Maho  All rights reserved.
 * BSD 2-Clause License (see mplapack_arithmetic_params.h for full text).
 */

#pragma once
#ifndef MPLAPACK_ARITHMETIC_PARAMS_DOUBLE_H
#define MPLAPACK_ARITHMETIC_PARAMS_DOUBLE_H

// Precondition: mplapack_arithmetic_params.h has been included.
// Precondition: <limits> is available.

#include <limits>

namespace mplapack {

// ---------------------------------------------------------------------------
// ArithmeticParams<double>
//
// Machine metadata is taken from std::numeric_limits<double>.
// For IEEE binary64:
//   digits = 53
//   emin   = -1021
//   emax   = 1024
// ---------------------------------------------------------------------------
template <> inline ArithmeticParams<double> get_arithmetic_params<double>() {
    ArithmeticParams<double> p;

    using limits = std::numeric_limits<double>;

    // Rlamch("E"): unit roundoff for round-to-nearest = 2^(-digits)
    const double eps = limits::epsilon() / 2.0;

    p.eps = eps;
    p.base = 2.0;
    p.prec = eps * p.base;
    p.t = static_cast<arithmetic_int>(limits::digits);
    p.rnd = 1.0; // IEEE 754 round-to-nearest
    p.emin = static_cast<arithmetic_int>(limits::min_exponent);
    p.rmin = limits::min();
    p.emax = static_cast<arithmetic_int>(limits::max_exponent);
    p.rmax = limits::max();

    // Rlamch("S"): safe minimum following netlib DLAMCH('S')
    p.sfmin = limits::min();
    const double small = 1.0 / limits::max();
    if (small >= p.sfmin) {
        p.sfmin = small * (1.0 + eps);
    }

    p.safmin = detail::compute_safmin<double>(p.emin, p.emax);
    p.safmax = detail::compute_safmax<double>(p.emin, p.emax);

    return p;
}

// ---------------------------------------------------------------------------
// BlueScalingParams<double>
//
// Reference exponents for IEEE binary64:
//   emin = -1021, emax = 1024, digits = 53
//   exp_tsml = ceil((-1021-1)/2)        = -511
//   exp_tbig = floor((1024-53+1)/2)     =  486
//   exp_ssml = -floor((-1021-53)/2)     =  537
//   exp_sbig = -ceil((1024+53-1)/2)     = -538
//
// Values are derived via make_blue_scaling_params(get_arithmetic_params<T>()).
// ---------------------------------------------------------------------------
template <> inline BlueScalingParams<double> get_blue_scaling_params<double>() {
    return make_blue_scaling_params(get_arithmetic_params<double>());
}

} // namespace mplapack
#endif // MPLAPACK_ARITHMETIC_PARAMS_DOUBLE_H

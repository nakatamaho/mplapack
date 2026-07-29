/*
 * Copyright (c) 2008-2025
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
 *
 */

// Mexponent: Fortran EXPONENT intrinsic for all MPLAPACK backends.
//
// Fortran model: x = f * 2^e, where 0.5 <= |f| < 1.  Returns e.
//
// C/C++ ilogb uses x = f * 2^e with 1 <= |f| < 2, so
//   Fortran EXPONENT(x) == std::ilogb(x) + 1.
// mpfr_get_exp() and mpf_get_d_2exp() already follow the 0.5 <= |f| < 1
// model, so they return the Fortran EXPONENT value without adjustment.

#include <mpblas.h>
#include <mplapack.h>

#include <cmath>

// -----------------------------------------------------------------------
// double
// -----------------------------------------------------------------------
#if defined MPLAPACK_BUILD_WITH_DOUBLE
INTEGER Mexponent(REAL const &x) {
    return std::ilogb(x) + 1;
}
#endif

// -----------------------------------------------------------------------
// binary80
// -----------------------------------------------------------------------
#if defined MPLAPACK_BUILD_WITH_BINARY80
#  if MPLAPACK_BINARY80_MATH == MPLAPACK_BINARY80_MATH_LDBL
// long double == binary80 (x87 on x86/x86-64).
INTEGER Mexponent(REAL const &x) {
    return std::ilogbl(x) + 1;
}
#  elif MPLAPACK_BINARY80_MATH == MPLAPACK_BINARY80_MATH_F64X
// _Float64x: cast to long double is exact when _Float64x == binary80.
INTEGER Mexponent(REAL const &x) {
    return std::ilogbl((long double)x) + 1;
}
#  else
#    error "Mexponent: unsupported MPLAPACK_BINARY80_MATH"
#  endif
#endif // MPLAPACK_BUILD_WITH_BINARY80

// -----------------------------------------------------------------------
// binary128
// -----------------------------------------------------------------------
#if defined MPLAPACK_BUILD_WITH_BINARY128
#  if MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_LDBL
// long double == binary128 (SPARC, MIPS, some AArch64).
INTEGER Mexponent(REAL const &x) {
    return std::ilogbl(x) + 1;
}
#  elif MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_F128
// _Float128 (ISO C TS 18661-3, GCC 7+).
// __builtin_ilogbf128 accepts _Float128 directly on GCC.
INTEGER Mexponent(REAL const &x) {
    return __builtin_ilogbf128((_Float128)x) + 1;
}
#  elif MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_QUADMATH
// _Float128 backed by libquadmath.  ilogbq expects __float128, but
// _Float128 and __float128 are the same underlying type on GCC, so the
// cast is a no-op.
#    include <quadmath.h>
INTEGER Mexponent(REAL const &x) {
    return ilogbq((__float128)x) + 1;
}
#  else
#    error "Mexponent: unsupported MPLAPACK_BINARY128_MATH"
#  endif
#endif // MPLAPACK_BUILD_WITH_BINARY128

// -----------------------------------------------------------------------
// DD (double-double)
// -----------------------------------------------------------------------
#if defined MPLAPACK_BUILD_WITH_DD
INTEGER Mexponent(REAL const &x) {
    // x.x[0] is the dominant component; lower limbs are bounded by
    // ulp(x.x[0]) and can never alter the binade.
    return std::ilogb(x.x[0]) + 1;
}
#endif // MPLAPACK_BUILD_WITH_DD

// -----------------------------------------------------------------------
// QD (quad-double)
// -----------------------------------------------------------------------
#if defined MPLAPACK_BUILD_WITH_QD
INTEGER Mexponent(REAL const &x) {
    // Same reasoning as DD.
    return std::ilogb(x.x[0]) + 1;
}
#endif // MPLAPACK_BUILD_WITH_QD

// -----------------------------------------------------------------------
// MPFR (mpfr::mpreal)
// -----------------------------------------------------------------------
#if defined MPLAPACK_BUILD_WITH_MPFR
INTEGER Mexponent(REAL const &x) {
    // mpfr_get_exp returns e s.t. 0.5 <= |x| * 2^(-e) < 1:
    // exactly Fortran EXPONENT semantics.  No +1 needed.
    return static_cast<INTEGER>(mpfr_get_exp((mpfr_ptr)const_cast<REAL &>(x)));
}
#endif // MPLAPACK_BUILD_WITH_MPFR

// -----------------------------------------------------------------------
// GMP (mpf_class)
// -----------------------------------------------------------------------
#if defined MPLAPACK_BUILD_WITH_GMP
INTEGER Mexponent(REAL const &x) {
    // mpf_get_d_2exp: sets exp s.t. x = d * 2^exp, 0.5 <= |d| < 1.
    // exp is the Fortran EXPONENT value directly.  No +1 needed.
    long exp;
    mpf_get_d_2exp(&exp, x.get_mpf_t());
    return static_cast<INTEGER>(exp);
}
#endif // MPLAPACK_BUILD_WITH_GMP

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

// Risinf: Inf predicate for all MPLAPACK backends.
//
// Companion to Risnan. Returns true if x is +Inf or -Inf.
// IEEE 754 isinf semantics are used for hardware types.
// For multi-component types (DD, QD), use the QD library predicates.
// For GMP (mpf_class), Inf is not representable; always returns false.

#include <mpblas.h>
#include <mplapack.h>

#include <cmath>

// -----------------------------------------------------------------------
// double
// -----------------------------------------------------------------------
#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
bool Risinf(REAL const &x) { return std::isinf(x); }
#endif

// -----------------------------------------------------------------------
// binary80
// -----------------------------------------------------------------------
#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
#if MPLAPACK_BINARY80_MATH == MPLAPACK_BINARY80_MATH_LDBL
bool Risinf(REAL const &x) { return std::isinf(x); }
#elif MPLAPACK_BINARY80_MATH == MPLAPACK_BINARY80_MATH_F64X
// _Float64x: cast to long double is exact when _Float64x == binary80.
bool Risinf(REAL const &x) { return std::isinf((long double)x); }
#else
#error "Risinf: unsupported MPLAPACK_BINARY80_MATH"
#endif
#endif // ___MPLAPACK_BUILD_WITH_BINARY80___

// -----------------------------------------------------------------------
// binary128
// -----------------------------------------------------------------------
#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
#if MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_LDBL
// long double == binary128 (SPARC, MIPS, some AArch64).
bool Risinf(REAL const &x) { return std::isinf(x); }
#elif MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_F128
// _Float128 (ISO C TS 18661-3, GCC 7+).
bool Risinf(REAL const &x) { return __builtin_isinf((_Float128)x); }
#elif MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_QUADMATH
#include <quadmath.h>
bool Risinf(REAL const &x) { return isinfq((__float128)x); }
#else
#error "Risinf: unsupported MPLAPACK_BINARY128_MATH"
#endif
#endif // ___MPLAPACK_BUILD_WITH_BINARY128___

// -----------------------------------------------------------------------
// DD (double-double)
// -----------------------------------------------------------------------
#if defined ___MPLAPACK_BUILD_WITH_DD___
bool Risinf(REAL const &x) {
    return x.isinf();
}
#endif // ___MPLAPACK_BUILD_WITH_DD___

// -----------------------------------------------------------------------
// QD (quad-double)
// -----------------------------------------------------------------------
#if defined ___MPLAPACK_BUILD_WITH_QD___
bool Risinf(REAL const &x) {
    return x.isinf();
}
#endif // ___MPLAPACK_BUILD_WITH_QD___

// -----------------------------------------------------------------------
// MPFR (mpfr::mpreal)
// -----------------------------------------------------------------------
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
bool Risinf(REAL const &x) {
    return mpfr_inf_p((mpfr_ptr)const_cast<REAL &>(x)) != 0;
}
#endif // ___MPLAPACK_BUILD_WITH_MPFR___

// -----------------------------------------------------------------------
// GMP (mpf_class)
// -----------------------------------------------------------------------
#if defined ___MPLAPACK_BUILD_WITH_GMP___
bool Risinf(REAL const &x) {
    // GMP mpf does not support NaN or Inf; all mpf values are finite.
    // Always return false.
    (void)x;
    return false;
}
#endif // ___MPLAPACK_BUILD_WITH_GMP___

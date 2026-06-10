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
 *
 */

// ---------------------------------------------------------------------------
// Rlamch.cpp -- thin adapter layer
//
// For MPFR, GMP, double, binary80, binary128: all logic lives in the
// mplapack_arithmetic_params_<type>.h headers; this file contains only the
// character-dispatch wrapper that forwards to get_arithmetic_params<REAL>().
//
// For QD and DD: the original implementations are retained verbatim until
// mplapack_arithmetic_params_qd.h / mplapack_arithmetic_params_dd.h are added.
// ---------------------------------------------------------------------------

#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_arithmetic_params.h>

#include <stdio.h>

#pragma STDC FP_CONTRACT OFF
#pragma STDC FENV_ACCESS ON

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("O0")
#pragma GCC optimize("no-tree-vectorize")
#endif

// =============================================================================
// Backends that use the new arithmetic_params layer.
// A single macro-guarded block covers all five ported types; each produces
// a Rlamch_<type>(cmach) function that maps to the appropriate struct field.
// =============================================================================

template <class T>
inline T rlamch_dispatch_impl(const char *cmach) {
    const auto p = mplapack::get_arithmetic_params<T>();
    if (Mlsame(cmach, "E")) return p.eps;
    if (Mlsame(cmach, "S")) return p.sfmin;
    if (Mlsame(cmach, "B")) return p.base;
    if (Mlsame(cmach, "P")) return p.prec;
    if (Mlsame(cmach, "N")) return mplapack::detail::to_rlamch_real<T>(p.t);
    if (Mlsame(cmach, "R")) return p.rnd;
    if (Mlsame(cmach, "M")) return mplapack::detail::to_rlamch_real<T>(p.emin);
    if (Mlsame(cmach, "U")) return p.rmin;
    if (Mlsame(cmach, "L")) return mplapack::detail::to_rlamch_real<T>(p.emax);
    if (Mlsame(cmach, "O")) return p.rmax;
    Mxerbla("Rlamch", 1);
    return mplapack::detail::to_rlamch_real<T>(0);
}

// ---------------------------------------------------------------------------
// MPFR
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_MPFR___)
REAL Rlamch_mpfr(const char *cmach) {
    return rlamch_dispatch_impl<mpreal>(cmach);
}
#endif

// ---------------------------------------------------------------------------
// GMP
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_GMP___)
REAL Rlamch_gmp(const char *cmach) {
    return rlamch_dispatch_impl<mpf_class>(cmach);
}
#endif

// ---------------------------------------------------------------------------
// double
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_DOUBLE___)
double Rlamch_double(const char *cmach) {
    return rlamch_dispatch_impl<double>(cmach);
}
#endif

// ---------------------------------------------------------------------------
// binary80
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_BINARY80___)
mplapack_binary80_t Rlamch_binary80(const char *cmach) {
    return rlamch_dispatch_impl<mplapack_binary80_t>(cmach);
}
#endif

// ---------------------------------------------------------------------------
// binary128
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_BINARY128___)
mplapack_binary128_t Rlamch_binary128(const char *cmach) {
    return rlamch_dispatch_impl<mplapack_binary128_t>(cmach);
}
#endif

// ---------------------------------------------------------------------------
// QD
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_QD___)
qd_real Rlamch_qd(const char *cmach) {
    return rlamch_dispatch_impl<qd_real>(cmach);
}
#endif

// ---------------------------------------------------------------------------
// DD
// ---------------------------------------------------------------------------
#if defined(___MPLAPACK_BUILD_WITH_DD___)
dd_real Rlamch_dd(const char *cmach) {
    return rlamch_dispatch_impl<dd_real>(cmach);
}
#endif

#undef MPLAPACK_RLAMCH_DISPATCH

// =============================================================================
// Rlamc3: always-present
// =============================================================================
REAL Rlamc3(REAL a, REAL b) { return a + b; }

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

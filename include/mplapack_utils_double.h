/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_utils_double.h,v 1.11 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MUTILS_DOUBLE_H_
#define _MUTILS_DOUBLE_H_

#if !(_XOPEN_SOURCE >= 600) && !(_ISOC99_SOURCE)
inline double log2(double x) { return log(x) / log(2.0); }
#endif

#if defined ___MPLAPACK_INTERNAL___
#include <mplapack_print_double.h>
#endif

inline double pow2(double a) {
    double mtmp;
    mtmp = a * a;
    return mtmp;
}

#include <type_traits>

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

// implementation of sign transfer function.
inline double sign(double a, double b) {
    double mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline double castREAL_double(mplapackint n) {
    double ret = n;
    return ret;
}

inline mplapackint castINTEGER_double(double a) {
    mplapackint i = a;
    return i;
}

inline long nint(double a) {
    long i;
    double tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (long)tmp;
    return i;
}

inline double cast2double(double a) { return a; }

inline std::complex<double> exp(std::complex<double> x) {
    double ex;
    double c;
    double s;
    std::complex<double> ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

inline double pi(double dummy) { return M_PI; }

static inline double cabs1(const std::complex<double> &z) { return abs(z.real()) + abs(z.imag()); }

#include <type_traits>

// NOTE:
// Do NOT 'using std::min/max' here.
// std::min/max have a 3-arg overload where the 3rd argument is a comparator,
// which hijacks Fortran-style min(a,b,c)/max(a,b,c) calls.

#ifndef MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED

// min/max for mplapackint (Fortran INTEGER)
inline mplapackint min(mplapackint a, mplapackint b) { return (a > b) ? b : a; }
inline mplapackint max(mplapackint a, mplapackint b) { return (a < b) ? b : a; }

// 3-arg overloads block std::min/max(a,b,comp) hijack.
inline mplapackint min(mplapackint a, mplapackint b, mplapackint c) {
    mplapackint r = min(a, b);
    return min(r, c);
}
inline mplapackint max(mplapackint a, mplapackint b, mplapackint c) {
    mplapackint r = max(a, b);
    return max(r, c);
}

// 4+ args: fold expression, mplapackint only.
template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mplapackint, std::decay_t<Args>> && ...)>> inline mplapackint min(mplapackint a, mplapackint b, mplapackint c, Args... rest) {
    mplapackint r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mplapackint, std::decay_t<Args>> && ...)>> inline mplapackint max(mplapackint a, mplapackint b, mplapackint c, Args... rest) {
    mplapackint r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_MPLAPACKINT_DEFINED

#include <type_traits>

#ifndef MPLAPACK_MINMAX_DOUBLE_DEFINED
#define MPLAPACK_MINMAX_DOUBLE_DEFINED

inline double min(double a, double b) { return (a > b) ? b : a; }
inline double max(double a, double b) { return (a < b) ? b : a; }

inline double min(double a, double b, double c) {
    double r = min(a, b);
    return min(r, c);
}
inline double max(double a, double b, double c) {
    double r = max(a, b);
    return max(r, c);
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<double, std::decay_t<Args>> && ...)>> inline double min(double a, double b, double c, Args... rest) {
    double r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<double, std::decay_t<Args>> && ...)>> inline double max(double a, double b, double c, Args... rest) {
    double r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_DOUBLE_DEFINED

#ifndef MPLAPACK_CHAR_UTILS_H
#define MPLAPACK_CHAR_UTILS_H

// Small helpers to build short option strings for ILAENV / IMlaenv calls.
//
// Typical usage:
//   const char *jbcmpz = CHAR2(job, compz);
//   nmin = iMlaenv(12, "Chseqr", jbcmpz, n, ilo, ihi, lwork);
//
// job, compz, side, howmny, sense, ... are almost always const char* pointing
// to single-character flags ("N", "V", "S", "E", etc.).

// 2-character helper -------------------------------------------------------
inline const char *CHAR2(char c1, char c2) {
    // Thread-local to avoid cross-call races.
    static thread_local char buf[3];
    buf[0] = c1;
    buf[1] = c2;
    buf[2] = '\0';
    return buf;
}

inline const char *CHAR2(const char *c1, const char *c2) {
    // Accept "N", "V", etc. as const char* and take their first characters.
    const char a = (c1 && *c1) ? *c1 : '\0';
    const char b = (c2 && *c2) ? *c2 : '\0';
    return CHAR2(a, b);
}

// 3-character helper -------------------------------------------------------
inline const char *CHAR3(char c1, char c2, char c3) {
    static thread_local char buf[4];
    buf[0] = c1;
    buf[1] = c2;
    buf[2] = c3;
    buf[3] = '\0';
    return buf;
}

inline const char *CHAR3(const char *c1, const char *c2, const char *c3) {
    const char a = (c1 && *c1) ? *c1 : '\0';
    const char b = (c2 && *c2) ? *c2 : '\0';
    const char c = (c3 && *c3) ? *c3 : '\0';
    return CHAR3(a, b, c);
}

#endif // MPLAPACK_CHAR_UTILS_H

#endif

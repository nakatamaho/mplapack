/*
 * Copyright (c) 2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
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

#ifndef _MPLAPACK_UTILS__FLOAT64X_H_
#define _MPLAPACK_UTILS__FLOAT64X_H_

#include "mplapack_config.h"

#if defined ___MPLAPACK_USE___FLOAT80___
// _Float64x is not available as a native type on this compiler; use __float80.
typedef __float80 _Float64x;
#endif

#include <math.h>

#if defined ___MPLAPACK__FLOAT64X_IS_LONGDOUBLE___

// Intel C/C++ compilers currently only supports long double as extended doule (aka FP80)
// and following workaround works
#if defined(__INTEL_COMPILER)
    typedef long double _Float64x;
#endif

#if defined ___MPLAPACK_INTERNAL___
#define FLOAT64X_FORMAT "%+25.21Le"
#define FLOAT64X_SHORT_FORMAT "%+20.16Le"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(_Float64x rtmp) { printf(FLOAT64X_FORMAT, rtmp); }
inline void printnum(std::complex<_Float64x> ctmp) { printf(FLOAT64X_FORMAT FLOAT64X_FORMAT "i", ctmp.real(), ctmp.imag()); }

inline void printnum_short(_Float64x rtmp) { printf(FLOAT64X_SHORT_FORMAT, rtmp); }
inline void printnum_short(std::complex<_Float64x> ctmp) { printf(FLOAT64X_SHORT_FORMAT FLOAT64X_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }

inline void sprintnum(char *buf, _Float64x rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_FORMAT, rtmp); }
inline void sprintnum(char *buf, std::complex<_Float64x> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_FORMAT FLOAT64X_FORMAT "i", ctmp.real(), ctmp.imag()); }

inline void sprintnum_short(char *buf, _Float64x rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_SHORT_FORMAT, rtmp); }
inline void sprintnum_short(char *buf, std::complex<_Float64x> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_SHORT_FORMAT FLOAT64X_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }
#endif

inline _Float64x pow2(const _Float64x &a) { return a * a; }
inline std::complex<_Float64x> pow2(const std::complex<_Float64x> &a) { return a * a; }

#ifdef __cplusplus
extern "C" {
#endif
#include <complex.h>
#include <complex>
#ifdef __cplusplus
}
#endif

// implementation of sign transfer function.
inline _Float64x sign(const _Float64x &a, const _Float64x &b) {
    _Float64x mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline double cast2double(_Float64x a) { return (double)a; }

inline long nint(_Float64x a) {
    long i;
    _Float64x tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (long)tmp;
    return i;
}

inline _Float64x castREAL__Float64x(mplapackint n) {
    _Float64x ret = n;
    return ret;
}

inline mplapackint castINTEGER__Float64x(_Float64x a) {
    mplapackint i = a;
    return i;
}

inline _Float64x pi(_Float64x dummy) {
#if defined __APPLE__ // __MATH_LONG_DOUBLE_CONSTANTS looks broken
    return 0xc.90fdaa22168c235p-2L;
#elif defined __MINGW32__
    return 0xc.90fdaa22168c235p-2L;
#else
    return M_PIl;
#endif
}

#endif

#if defined ___MPLAPACK_HAVE_NATIVE__FLOAT64X___

#include <float.h>

#define STR2(x) #x
#define STR(x) STR2(x)

#if defined(__FLT64X_MANT_DIG__) && defined(__FLT64X_MAX_EXP__) && defined(LDBL_MANT_DIG) && defined(LDBL_MAX_EXP)

// We can compare the formats of _Float64x and long double.
#if (__FLT64X_MANT_DIG__ == LDBL_MANT_DIG) && (__FLT64X_MAX_EXP__ == LDBL_MAX_EXP)
// Formats match -> OK
#define ___MPLAPACK__FLOAT64X_IS_LONGDOUBLE___ 1
#else
// Formats differ -> print all four numbers
#error "_Float64x != long double: __FLT64X_MANT_DIG__=" \
           STR(__FLT64X_MANT_DIG__)                       \
           " LDBL_MANT_DIG=" STR(LDBL_MANT_DIG)          \
           " __FLT64X_MAX_EXP__=" STR(__FLT64X_MAX_EXP__) \
           " LDBL_MAX_EXP=" STR(LDBL_MAX_EXP)
#endif
#else
#error "Cannot compare _Float64x and long double: missing macros: " \
         "_FLT64X_MANT_DIG__=" STR(__FLT64X_MANT_DIG__) " " \
         "_FLT64X_MAX_EXP__=" STR(__FLT64X_MAX_EXP__) " "  \
         "LDBL_MANT_DIG=" STR(LDBL_MANT_DIG) " "            \
         "LDBL_MAX_EXP=" STR(LDBL_MAX_EXP)
#endif

// Intel C/C++ compilers currently only supports long double as extended doule (aka FP80)
// and following workaround works
#if defined(__INTEL_COMPILER)
typedef long double _Float64x;
#endif

#if defined ___MPLAPACK_INTERNAL___
#define FLOAT64X_FORMAT "%+25.21Le"
#define FLOAT64X_SHORT_FORMAT "%+20.16Le"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(_Float64x rtmp) { printf(FLOAT64X_FORMAT, (long double)rtmp); }
inline void printnum(std::complex<_Float64x> ctmp) { printf(FLOAT64X_FORMAT FLOAT64X_FORMAT "i", (long double)ctmp.real(), (long double)ctmp.imag()); }

inline void printnum_short(_Float64x rtmp) { printf(FLOAT64X_SHORT_FORMAT, (long double)rtmp); }
inline void printnum_short(std::complex<_Float64x> ctmp) { printf(FLOAT64X_SHORT_FORMAT FLOAT64X_SHORT_FORMAT "i", (long double)ctmp.real(), (long double)ctmp.imag()); }

inline void sprintnum(char *buf, _Float64x rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_FORMAT, (long double)rtmp); }
inline void sprintnum(char *buf, std::complex<_Float64x> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_FORMAT FLOAT64X_FORMAT "i", (long double)ctmp.real(), (long double)ctmp.imag()); }

inline void sprintnum_short(char *buf, _Float64x rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_SHORT_FORMAT, (long double)rtmp); }
inline void sprintnum_short(char *buf, std::complex<_Float64x> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT64X_SHORT_FORMAT FLOAT64X_SHORT_FORMAT "i", (long double)ctmp.real(), (long double)ctmp.imag()); }
#endif

inline _Float64x pow2(const _Float64x &a) { return a * a; }

#ifdef __cplusplus
extern "C" {
#endif
#include <complex.h>
#include <complex>
#ifdef __cplusplus
}
#endif

// implementation of sign transfer function.
inline _Float64x sign(const _Float64x &a, const _Float64x &b) {
    _Float64x mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline double cast2double(_Float64x a) { return (double)a; }

inline long nint(_Float64x a) {
    long i;
    _Float64x tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (long)tmp;
    return i;
}

inline _Float64x castREAL__Float64x(mplapackint n) {
    _Float64x ret = n;
    return ret;
}

inline mplapackint castINTEGER__Float64x(_Float64x a) {
    mplapackint i = a;
    return i;
}

inline _Float64x pi(_Float64x dummy) {
#if defined __APPLE__ // __MATH_LONG_DOUBLE_CONSTANTS looks broken
    return 0xc.90fdaa22168c235p-2L;
#elif defined __MINGW32__
    return 0xc.90fdaa22168c235p-2L;
#else
    return M_PIl;
#endif
}

#endif

#if defined ___MPLAPACK__FLOAT64X_BROKEN___
inline _Float64x abs(const std::complex<_Float64x> &a) { return sqrtl(a.real() * a.real() + a.imag() * a.imag()); }

inline std::complex<_Float64x> sqrt(const std::complex<_Float64x> a) {
    _Float64x _Complex b, tmp;
    std::complex<_Float64x> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csqrtl(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float64x> sin(const std::complex<_Float64x> a) {
    _Float64x _Complex b, tmp;
    std::complex<_Float64x> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csinl(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float64x> cos(const std::complex<_Float64x> a) {
    _Float64x _Complex b, tmp;
    std::complex<_Float64x> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = ccosl(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float64x> exp(const std::complex<_Float64x> &a) {
    _Float64x _Complex b, tmp;
    std::complex<_Float64x> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = cexpl(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float64x> log(const std::complex<_Float64x> &a) {
    _Float64x _Complex b, tmp;
    std::complex<_Float64x> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = clogl(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}
#endif //___MPLAPACK__FLOAT64X_BROKEN___

static inline _Float64x cabs1(const std::complex<_Float64x> &z) { return abs(z.real()) + abs(z.imag()); }

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

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

#include <type_traits>

#ifndef MPLAPACK_MINMAX_FLOAT64X_DEFINED
#define MPLAPACK_MINMAX_FLOAT64X_DEFINED

inline _Float64x min(_Float64x a, _Float64x b) { return (a > b) ? b : a; }
inline _Float64x max(_Float64x a, _Float64x b) { return (a < b) ? b : a; }

inline _Float64x min(_Float64x a, _Float64x b, _Float64x c) {
    _Float64x r = min(a, b);
    return min(r, c);
}
inline _Float64x max(_Float64x a, _Float64x b, _Float64x c) {
    _Float64x r = max(a, b);
    return max(r, c);
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<_Float64x, std::decay_t<Args>> && ...)>> inline _Float64x min(_Float64x a, _Float64x b, _Float64x c, Args... rest) {
    _Float64x r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<_Float64x, std::decay_t<Args>> && ...)>> inline _Float64x max(_Float64x a, _Float64x b, _Float64x c, Args... rest) {
    _Float64x r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_FLOAT64X_DEFINED

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

// Integer ceil for _Float64x.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_FLOAT64X_DEFINED
#define MPLAPACK_ICEIL_FLOAT64X_DEFINED
inline mplapackint iceil(_Float64x x) {
    mplapackint t = static_cast<mplapackint>(x); // trunc toward zero
    if (x > static_cast<_Float64x>(t)) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_FLOAT64X_DEFINED

#endif

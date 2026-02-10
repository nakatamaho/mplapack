/*
 * Copyright (c) 2012-2026
 *	Nakata, Maho
 * 	All rights reserved.
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

#ifndef _MUTILS__BINARY128_H_
#define _MUTILS__BINARY128_H_

#include "mplapack_config.h"

#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#include <quadmath.h>
#endif // MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH

#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
#if defined ___MPLAPACK_INTERNAL___
#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif
#include <string.h>

#define BINARY128_FORMAT "%+-#*.40Qe"
#define BINARY128_SHORT_FORMAT "%+-#*.16Qe"

inline void printnum(__float128 rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__];
    int n = quadmath_snprintf(buf, sizeof buf, BINARY128_FORMAT, width, rtmp);
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    return;
}

inline void printnum(std::complex<__float128> rtmp) {
    int width = 42, n;
    char buf[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    n = quadmath_snprintf(buf, sizeof buf, BINARY128_FORMAT, width, rtmp.real());
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    n = quadmath_snprintf(buf2, sizeof buf, BINARY128_FORMAT, width, rtmp.imag());
    if ((size_t)n < sizeof buf2)
        printf("%s", buf2);
    printf("i");
    return;
}

inline void sprintnum(char *buf, __float128 rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, BINARY128_FORMAT, width, rtmp);
    return;
}

inline void sprintnum(char *buf, std::complex<__float128> rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, BINARY128_FORMAT BINARY128_FORMAT, width, rtmp.real(), rtmp.imag());
    return;
}

inline void printnum_short(__float128 rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__];
    int n = quadmath_snprintf(buf, sizeof buf, BINARY128_SHORT_FORMAT, width, rtmp);
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    return;
}

inline void printnum_short(std::complex<__float128> rtmp) {
    int width = 42, n;
    char buf[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    n = quadmath_snprintf(buf, sizeof buf, BINARY128_SHORT_FORMAT, width, rtmp.real());
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    n = quadmath_snprintf(buf2, sizeof buf, BINARY128_SHORT_FORMAT, width, rtmp.imag());
    if ((size_t)n < sizeof buf2)
        printf("%s", buf2);
    printf("i");
    return;
}

inline void sprintnum_short(char *buf, __float128 rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, BINARY128_SHORT_FORMAT, width, rtmp);
    return;
}

inline void sprintnum_short(char *buf, std::complex<__float128> rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, BINARY128_SHORT_FORMAT BINARY128_SHORT_FORMAT, width, rtmp.real(), rtmp.imag());
    return;
}
#endif //___MPLAPACK_INTERNAL___
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
#if defined ___MPLAPACK_INTERNAL___
#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif
#include <string.h>

// Full precision format (40 digits)
#define FLOAT128_FORMAT "%.40e"
// Short format (16 digits)
#define FLOAT128_SHORT_FORMAT "%.16e"

// printnum - full precision output to stdout
inline void printnum(_Float128 rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__];
    strfromf128(buf, sizeof buf, FLOAT128_FORMAT, rtmp);
    printf("%s", buf);
    return;
}

inline void printnum(std::complex<_Float128> rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    strfromf128(buf, sizeof buf, FLOAT128_FORMAT, rtmp.real());
    printf("%s", buf);
    strfromf128(buf2, sizeof buf, FLOAT128_FORMAT, rtmp.imag());
    printf("%s", buf2);
    printf("i");
    return;
}

// sprintnum - full precision output to buffer
inline void sprintnum(char *buf, _Float128 rtmp) {
    int width = 42;
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp);
    return;
}

inline void sprintnum(char *buf, std::complex<_Float128> rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    buf[0] = '\0';
    if (rtmp.real() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf1, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp.real());
    strcat(buf, buf1);
    if (rtmp.imag() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf2, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp.imag());
    strcat(buf, buf2);
    strcat(buf, "i");
    return;
}

// printnum_short - short precision output to stdout
inline void printnum_short(_Float128 rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__];
    strfromf128(buf, sizeof buf, FLOAT128_SHORT_FORMAT, rtmp);
    printf("%s", buf);
    return;
}

inline void printnum_short(std::complex<_Float128> rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    strfromf128(buf, sizeof buf, FLOAT128_SHORT_FORMAT, rtmp.real());
    printf("%s", buf);
    strfromf128(buf2, sizeof buf, FLOAT128_SHORT_FORMAT, rtmp.imag());
    printf("%s", buf2);
    printf("i");
    return;
}

// sprintnum_short - short precision output to buffer
inline void sprintnum_short(char *buf, _Float128 rtmp) {
    int width = 42;
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp);
    return;
}

inline void sprintnum_short(char *buf, std::complex<_Float128> rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    buf[0] = '\0';
    if (rtmp.real() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf1, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp.real());
    strcat(buf, buf1);
    if (rtmp.imag() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf2, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp.imag());
    strcat(buf, buf2);
    strcat(buf, "i");
    return;
}

#endif // ___MPLAPACK_INTERNAL___
#else
#error
#endif // MPLAPACK_BINARY128_IO

#if MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_QUADMATH
#include <quadmath.h>

#ifndef MPLAPACK_HAVE_STD_ABS_FLOAT128
#define MPLAPACK_HAVE_STD_ABS_FLOAT128 0
#endif

inline __float128 pow(const __float128 &a, const __float128 &b) { return powq(a, b); }
inline __float128 pow(const long &a, const long &b) { return powq((__float128)a, (__float128)b); }
inline __float128 pow(const int &a, const long &b) { return powq((__float128)a, (__float128)b); }
inline __float128 pow(const __float128 &a, const long &b) { return powq(a, (__float128)b); }
inline __float128 sqrt(const __float128 &a) { return sqrtq(a); }

#if !MPLAPACK_HAVE_STD_ABS_FLOAT128
// Define a fallback abs for __float128 only when std::abs(__float128) is missing.
inline __float128 abs(const __float128 &a) { return absq(a); }
#endif

inline __float128 sin(__float128 a) { return sinq(a); }
inline __float128 sinh(__float128 a) { return sinhq(a); }
inline __float128 cos(__float128 a) { return cosq(a); }
inline __float128 cosh(__float128 a) { return coshq(a); }
inline __float128 atan2(__float128 a, __float128 b) { return atan2q(a, b); }
inline __float128 exp(const __float128 &a) { return expq(a); }
inline __float128 log(const __float128 &a) { return logq(a); }
inline __float128 log10(const __float128 &a) { return log10q(a); }
inline __float128 log2(const __float128 &a) { return logq(a) / logq(2.0); }
inline __float128 ldexp(const __float128 &a, int exp) { return ldexpq(a, exp); }
inline __float128 nextafter(const __float128 &a, const __float128 &b) { return nextafterq(a, b); }
inline __float128 abs(const std::complex<__float128> &a) {
    __float128 _Complex b;
    __float128 c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    c = cabsq(b);
    return c;
}
inline std::complex<__float128> sqrt(const std::complex<__float128> a) {
    __float128 _Complex b, tmp;
    std::complex<__float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csqrtq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<__float128> sin(const std::complex<__float128> a) {
    __float128 _Complex b, tmp;
    std::complex<__float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csinq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<__float128> cos(const std::complex<__float128> a) {
    __float128 _Complex b, tmp;
    std::complex<__float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = ccosq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<__float128> exp(const std::complex<__float128> &a) {
    __float128 _Complex b, tmp;
    std::complex<__float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = cexpq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<__float128> log(const std::complex<__float128> &a) {
    __float128 _Complex b, tmp;
    std::complex<__float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = clogq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline long nint(mplapack_binary128_t a) {
    long i;
    mplapack_binary128_t tmp;
    a = a + 0.5;
    tmp = floorq(a);
    i = (long)tmp;
    return i;
}
inline __float128 pi(__float128 dummy) { return M_PIq; }

#elif MPLAPACK_BINARY128_MATH == MPLAPACK_BINARY128_MATH_F128
// Request IEC 60559 / TS 18661 (C23) extended floating and math/complex prototypes from glibc.
// Must be defined before including <math.h>/<complex.h> through <cmath>/<complex>.
#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif

#include <complex>
#include <complex.h>

inline _Float128 pow(const _Float128 &a, const _Float128 &b) { return powf128(a, b); }
inline _Float128 pow(const long &a, const long &b) { return powf128((_Float128)a, (_Float128)b); }
inline _Float128 pow(const int &a, const long &b) { return powf128((_Float128)a, (_Float128)b); }
inline _Float128 pow(const _Float128 &a, const long &b) { return powf128(a, (_Float128)b); }
inline _Float128 sqrt(const _Float128 &a) { return sqrtf128(a); }
inline _Float128 abs(const _Float128 &a) { return fabsf128(a); }

inline _Float128 sin(_Float128 a) { return sinf128(a); }
inline _Float128 sinh(_Float128 a) { return sinhf128(a); }
inline _Float128 cos(_Float128 a) { return cosf128(a); }
inline _Float128 cosh(_Float128 a) { return coshf128(a); }

inline _Float128 atan2(_Float128 a, _Float128 b) { return atan2f128(a, b); }

inline _Float128 exp(const _Float128 &a) { return expf128(a); }
inline _Float128 log(const _Float128 &a) { return logf128(a); }
inline _Float128 log10(const _Float128 &a) { return log10f128(a); }
inline _Float128 log2(const _Float128 &a) { return logf128(a) / logf128(2.0); }
inline _Float128 ceil(_Float128 a) { return ceilf128(a); }
inline _Float128 ldexp(const _Float128 &a, int exp) { return ldexpf128(a, exp); }

inline _Float128 abs(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    _Float128 c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    c = cabsf128(b);
    return c;
}

inline std::complex<_Float128> sqrt(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csqrtf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> sin(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csinf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> cos(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = ccosf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> exp(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = cexpf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> log(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = clogf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline mplapackint nint(_Float128 a) {
    mplapackint i;
    _Float128 tmp;
    a = a + 0.5;
    tmp = floorf128(a);
    i = (mplapackint)tmp;
    return i;
}

inline _Float128 pi(_Float128 dummy) { return M_PIf128; }
#endif // MPLAPACK_BINARY128_MATH

// support functions for mplapack, general ones
#include <cmath>

// implementation of sign transfer function.
inline mplapack_binary128_t sign(mplapack_binary128_t a, mplapack_binary128_t b) {
    mplapack_binary128_t mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline mplapack_binary128_t castREAL_binary128(mplapackint n) {
    mplapack_binary128_t ret = n;
    return ret;
}

inline mplapackint castINTEGER_binary128(mplapack_binary128_t a) {
    mplapackint i = a;
    return i;
}

inline double cast2double(mplapack_binary128_t a) { return (double)a; }

static inline mplapack_binary128_t cabs1(const std::complex<mplapack_binary128_t> &z) { return abs(z.real()) + abs(z.imag()); }

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

#ifndef MPLAPACK_MINMAX_BINARY128_DEFINED
#define MPLAPACK_MINMAX_BINARY128_DEFINED

inline mplapack_binary128_t min(mplapack_binary128_t a, mplapack_binary128_t b) { return (a > b) ? b : a; }
inline mplapack_binary128_t max(mplapack_binary128_t a, mplapack_binary128_t b) { return (a < b) ? b : a; }

inline mplapack_binary128_t min(mplapack_binary128_t a, mplapack_binary128_t b, mplapack_binary128_t c) {
    mplapack_binary128_t r = min(a, b);
    return min(r, c);
}
inline mplapack_binary128_t max(mplapack_binary128_t a, mplapack_binary128_t b, mplapack_binary128_t c) {
    mplapack_binary128_t r = max(a, b);
    return max(r, c);
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mplapack_binary128_t, std::decay_t<Args>> && ...)>> inline mplapack_binary128_t min(mplapack_binary128_t a, mplapack_binary128_t b, mplapack_binary128_t c, Args... rest) {
    mplapack_binary128_t r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mplapack_binary128_t, std::decay_t<Args>> && ...)>> inline mplapack_binary128_t max(mplapack_binary128_t a, mplapack_binary128_t b, mplapack_binary128_t c, Args... rest) {
    mplapack_binary128_t r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_BINARY128_DEFINED

#include <type_traits>

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

// Square and quartic for REAL (computational kernels)
#ifndef MPLAPACK_POW2_POW4_BINARY128_DEFINED
#define MPLAPACK_POW2_POW4_BINARY128_DEFINED
inline mplapack_binary128_t pow2(const mplapack_binary128_t &a) { return a * a; }
inline mplapack_binary128_t pow4(const mplapack_binary128_t &a) {
    mplapack_binary128_t tmp = a * a;
    return tmp * tmp;
}
inline std::complex<mplapack_binary128_t> pow2(const std::complex<mplapack_binary128_t> &a) { return a * a; }
inline std::complex<mplapack_binary128_t> pow4(const std::complex<mplapack_binary128_t> &a) {
    std::complex<mplapack_binary128_t> tmp = a * a;
    return tmp * tmp;
}
#endif // MPLAPACK_POW2_POW4_BINARY128_DEFINED

#ifndef MPLAPACK_CHAR_UTILS_H
#define MPLAPACK_CHAR_UTILS_H

// Small helpers to build short option strings for ILAENV / IMlaenv calls.
//
// Typical usage:
//   mnthr = iMlaenv(6, "Cgesvd", CHAR2(jobu, jobvt), m, n, 0, 0);
//
// jobu/jobvt/... are often const char* pointing to single-character flags
// ("N", "V", "S", "E", etc.).

struct charbuf2 {
    char s[3];
    constexpr charbuf2(char a, char b) : s{a, b, '\0'} {}
    constexpr operator const char *() const { return s; }
};

struct charbuf3 {
    char s[4];
    constexpr charbuf3(char a, char b, char c) : s{a, b, c, '\0'} {}
    constexpr operator const char *() const { return s; }
};

constexpr charbuf2 CHAR2(char a, char b) { return charbuf2(a, b); }
constexpr charbuf3 CHAR3(char a, char b, char c) { return charbuf3(a, b, c); }

// Extract first character from a 1-char C string (e.g. "N").
// If p is null or empty, returns '\0' to fail loudly downstream.
constexpr char first_char(const char *p) { return (p && p[0] != '\0') ? p[0] : '\0'; }

// Overloads for the common MPLAPACK/LAPACK style: const char* flags.
constexpr charbuf2 CHAR2(const char *a, const char *b) { return charbuf2(first_char(a), first_char(b)); }
constexpr charbuf3 CHAR3(const char *a, const char *b, const char *c) { return charbuf3(first_char(a), first_char(b), first_char(c)); }

#endif // MPLAPACK_CHAR_UTILS_H

// Integer ceil for mplapack_binary128_t.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_BINARY128_DEFINED
#define MPLAPACK_ICEIL_BINARY128_DEFINED
inline mplapackint iceil(mplapack_binary128_t x) {
    mplapackint t = static_cast<mplapackint>(x); // trunc toward zero
    if (x > static_cast<mplapack_binary128_t>(t)) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_BINARY128_DEFINED

#ifndef MPLAPACK_MOD_UTILS_H
#define MPLAPACK_MOD_UTILS_H
inline mplapackint mod(mplapackint a, mplapackint b) { return a % b; }
#endif // MPLAPACK_MOD_UTILS_H

#endif

/*
 * Copyright (c) 2021-2026
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

#ifndef _MPLAPACK_UTILS_BINARY80_H_
#define _MPLAPACK_UTILS_BINARY80_H_

#include "mplapack_config.h"

#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80

#endif // MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80

#if MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL
#if defined ___MPLAPACK_INTERNAL___
#define BINARY80_FORMAT "%+25.21Le"
#define BINARY80_SHORT_FORMAT "%+20.16Le"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(long double rtmp) { printf(BINARY80_FORMAT, rtmp); }
inline void printnum(std::complex<long double> ctmp) { printf(BINARY80_FORMAT BINARY80_FORMAT "i", ctmp.real(), ctmp.imag()); }

inline void printnum_short(long double rtmp) { printf(BINARY80_SHORT_FORMAT, rtmp); }
inline void printnum_short(std::complex<long double> ctmp) { printf(BINARY80_SHORT_FORMAT BINARY80_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }

inline void sprintnum(char *buf, long double rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, BINARY80_FORMAT, rtmp); }
inline void sprintnum(char *buf, std::complex<long double> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, BINARY80_FORMAT BINARY80_FORMAT "i", ctmp.real(), ctmp.imag()); }

inline void sprintnum_short(char *buf, long double rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, BINARY80_SHORT_FORMAT, rtmp); }
inline void sprintnum_short(char *buf, std::complex<long double> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, BINARY80_SHORT_FORMAT BINARY80_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }
#endif
#endif

#if MPLAPACK_BINARY80_MATH == MPLAPACK_BINARY80_MATH_LDBL

inline long double pow(const long double &a, const long double &b) { return powl(a, b); }
inline long double pow(const long &a, const long &b) { return powl((long double)a, (long double)b); }
inline long double pow(const int &a, const long &b) { return powl((long double)a, (long double)b); }
inline long double pow(const long double &a, const long &b) { return powl(a, (long double)b); }

inline long double sqrt(const long double &a) { return sqrtl(a); }
inline long double sin(long double a) { return sinl(a); }
inline long double sinh(long double a) { return sinhl(a); }
inline long double cos(long double a) { return cosl(a); }
inline long double cosh(long double a) { return coshl(a); }
inline long double atan2(long double a, long double b) { return atan2l(a, b); }
inline long double exp(const long double &a) { return expl(a); }
inline long double log(const long double &a) { return logl(a); }
inline long double log10(const long double &a) { return log10l(a); }
inline long double log2(const long double &a) { return logl(a) / logl(2.0L); }
inline long double ldexp(const long double &a, int exp) { return ldexpl(a, exp); }
inline long double nextafter(const long double &a, const long double &b) { return nextafterl(a, b); }

inline long double pow2(const long double &a) { return a * a; }
inline long double pow4(const long double &a) { return a * a * a * a; }
inline std::complex<long double> pow2(const std::complex<long double> &a) { return a * a; }
inline std::complex<long double> pow4(const std::complex<long double> &a) { return a * a * a * a; }

// implementation of sign transfer function.
inline long double sign(const long double &a, const long double &b) {
    long double mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline double cast2double(long double a) { return (double)a; }

inline long nint(long double a) {
    long i;
    long double tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (long)tmp;
    return i;
}

inline long double castREAL_binary80(mplapackint n) {
    long double ret = n;
    return ret;
}

inline mplapackint castINTEGER_binary80(long double a) {
    mplapackint i = a;
    return i;
}

inline long double pi(long double dummy) {
#if defined __APPLE__ // __MATH_LONG_DOUBLE_CONSTANTS looks broken
    return 0xc.90fdaa22168c235p-2L;
#elif defined __MINGW32__
    return 0xc.90fdaa22168c235p-2L;
#else
    return M_PIl;
#endif
}

#endif

static inline mplapack_binary80_t cabs1(const std::complex<mplapack_binary80_t> &z) { return abs(z.real()) + abs(z.imag()); }

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

#ifndef MPLAPACK_MINMAX_BINARY80_DEFINED
#define MPLAPACK_MINMAX_BINARY80_DEFINED

inline mplapack_binary80_t min(mplapack_binary80_t a, mplapack_binary80_t b) { return (a > b) ? b : a; }
inline mplapack_binary80_t max(mplapack_binary80_t a, mplapack_binary80_t b) { return (a < b) ? b : a; }

inline mplapack_binary80_t min(mplapack_binary80_t a, mplapack_binary80_t b, mplapack_binary80_t c) {
    mplapack_binary80_t r = min(a, b);
    return min(r, c);
}
inline mplapack_binary80_t max(mplapack_binary80_t a, mplapack_binary80_t b, mplapack_binary80_t c) {
    mplapack_binary80_t r = max(a, b);
    return max(r, c);
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mplapack_binary80_t, std::decay_t<Args>> && ...)>> inline mplapack_binary80_t min(mplapack_binary80_t a, mplapack_binary80_t b, mplapack_binary80_t c, Args... rest) {
    mplapack_binary80_t r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mplapack_binary80_t, std::decay_t<Args>> && ...)>> inline mplapack_binary80_t max(mplapack_binary80_t a, mplapack_binary80_t b, mplapack_binary80_t c, Args... rest) {
    mplapack_binary80_t r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_BINARY80_DEFINED

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

// Integer ceil for mplapack_binary80_t.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_BINARY80_DEFINED
#define MPLAPACK_ICEIL_BINARY80_DEFINED
inline mplapackint iceil(mplapack_binary80_t x) {
    mplapackint t = static_cast<mplapackint>(x); // trunc toward zero
    if (x > static_cast<mplapack_binary80_t>(t)) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_BINARY80_DEFINED

#ifndef MPLAPACK_MOD_UTILS_H
#define MPLAPACK_MOD_UTILS_H
inline mplapackint mod(mplapackint a, mplapackint b) { return a % b; }
#endif // MPLAPACK_MOD_UTILS_H

#endif

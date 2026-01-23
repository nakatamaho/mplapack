/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_utils_mpfr.h,v 1.9 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MUTILS_MPFR_H_
#define _MUTILS_MPFR_H_

#include "mpcomplex.h"
#include "mpreal.h"

using namespace mpfr;

#if defined ___MPLAPACK_INTERNAL___
#include <mplapack_print_double.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(mpreal rtmp) {
    mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp));
    return;
}

inline void printnum_short(mpreal rtmp) {
    mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void printnum(mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_FORMAT MPFR_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

inline void printnum_short(mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

inline void sprintnum(char *buf, mpreal rtmp) {
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum_short(char *buf, mpreal rtmp) {
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_SHORT_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum(char *buf, mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_FORMAT MPFR_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
inline void sprintnum_short(char *buf, mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
#endif

inline mpreal pow2(mpreal a) {
    mpreal mtmp = a * a;
    return mtmp;
}

inline mpcomplex pow2(mpcomplex a) {
    mpcomplex mtmp = a * a;
    return mtmp;
}

inline mpreal pow4(mpreal a) {
    mpreal mtmp = a * a * a *a;
    return mtmp;
}

inline mpcomplex pow4(mpcomplex a) {
    mpcomplex mtmp = a * a * a * a;
    return mtmp;
}

#include <type_traits>

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

// implementation of sign transfer function.
inline mpreal sign(mpreal a, mpreal b) {
    mpreal mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline mplapackint nint(mpreal a) {
    mplapackint i;
    mpreal tmp;
    a = a + 0.5;
    tmp = floorl(a);
    i = tmp; // cast to long
    return i;
}

inline mplapackint castINTEGER_mpfr(mpreal a) {
    mplapackint i;
    i = a;
    return i;
}

inline mpreal castREAL_mpfr(mplapackint a) {
    mpreal i = a;
    return i;
}

inline mpreal pi(mpreal dummy) {
    mpreal _PI;
    _PI = const_pi(mpfr::mpreal::default_prec);
    return _PI;
}

static inline mpreal cabs1(const mpcomplex &z) { return abs(z.real()) + abs(z.imag()); }

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

#ifndef MPLAPACK_MINMAX_MPREAL_VARIADIC_DEFINED
#define MPLAPACK_MINMAX_MPREAL_VARIADIC_DEFINED

// 3-arg overloads: blocks std::min/max(a,b,comp) hijack.
inline mpreal min(const mpreal &a, const mpreal &b, const mpreal &c) {
    mpreal r = (b < a) ? b : a;
    return (c < r) ? c : r;
}
inline mpreal max(const mpreal &a, const mpreal &b, const mpreal &c) {
    mpreal r = (a < b) ? b : a;
    return (r < c) ? c : r;
}

// 4+ args: fold expression, mpreal only.
template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mpreal, std::decay_t<Args>> && ...)>> inline mpreal min(const mpreal &a, const mpreal &b, const mpreal &c, const Args &...rest) {
    mpreal r = min(a, b, c);
    ((r = (rest < r) ? rest : r), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mpreal, std::decay_t<Args>> && ...)>> inline mpreal max(const mpreal &a, const mpreal &b, const mpreal &c, const Args &...rest) {
    mpreal r = max(a, b, c);
    ((r = (r < rest) ? rest : r), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_MPREAL_VARIADIC_DEFINED

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

// Integer ceil for MPFR mpreal.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_MPREAL_DEFINED
#define MPLAPACK_ICEIL_MPREAL_DEFINED
inline mplapackint iceil(const mpreal &x) {
    // mpreal -> integer cast truncates toward zero.
    mplapackint t = static_cast<mplapackint>(x);
    if (x > mpreal(t)) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_MPREAL_DEFINED

#endif

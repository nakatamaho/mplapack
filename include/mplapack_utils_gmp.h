/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_utils_gmp.h,v 1.10 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MUTILS_GMP_H_
#define _MUTILS_GMP_H_

#include "mpc_class.h"

#if defined ___MPLAPACK_INTERNAL___
#define GMP_FORMAT "%+68.64Fe"
#define GMP_SHORT_FORMAT "%+20.16Fe"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(mpf_class rtmp) {
    gmp_printf(GMP_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void printnum_short(mpf_class rtmp) {
    gmp_printf(GMP_SHORT_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void printnum(mpc_class ctmp) {
    gmp_printf(GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void printnum_short(mpc_class ctmp) {
    gmp_printf(GMP_SHORT_FORMAT GMP_SHORT_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void sprintnum(char *buf, mpf_class rtmp) {
    gmp_snprintf(buf, __MPLAPACK_BUFLEN__, GMP_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void sprintnum_short(char *buf, mpf_class rtmp) {
    gmp_snprintf(buf, __MPLAPACK_BUFLEN__, GMP_SHORT_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void sprintnum(char *buf, mpc_class ctmp) {
    gmp_snprintf(buf, __MPLAPACK_BUFLEN__, GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void sprintnum_short(char *buf, mpc_class ctmp) {
    gmp_snprintf(buf, __MPLAPACK_BUFLEN__, GMP_SHORT_FORMAT GMP_SHORT_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}
#endif

inline mpf_class pow2(mpf_class a) {
    mpf_class mtmp = a * a;
    return mtmp;
}

inline mpc_class pow2(mpc_class a) {
    mpc_class mtmp = a * a;
    return mtmp;
}

inline mpf_class pow4(mpf_class a) {
    mpf_class mtmp = a * a * a * a;
    return mtmp;
}

inline mpc_class pow4(mpc_class a) {
    mpc_class mtmp = a * a * a * a;
    return mtmp;
}

#include <type_traits>

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

inline mpf_class sign(mpf_class a, mpf_class b) {
    mpf_class mtmp;
    mpf_abs(mtmp.get_mpf_t(), a.get_mpf_t());
    if (b != 0.0) {
        mtmp = mpf_sgn(b.get_mpf_t()) * mtmp;
    }
    return mtmp;
}

inline mpf_class castREAL_gmp(mplapackint n) {
    mpf_class a(n);
    return a;
}

inline mplapackint castINTEGER_gmp(mpf_class a) {
    mplapackint i;
    i = mpf_get_si(a.get_mpf_t());
    return i;
}

inline mplapackint nint(mpf_class a) {
    mplapackint i;
    mpf_class tmp;
    a = a + 0.5;
    mpf_floor(tmp.get_mpf_t(), a.get_mpf_t());
    i = mpf_get_si(tmp.get_mpf_t());
    return i;
}

inline double cast2double(mpf_class a) { return a.get_d(); }

// every transcendental function and constant for GMP is in double precision.
inline mpf_class atan2(mpf_class a, mpf_class b) {
    double dtemp1, dtemp2;
    mpf_class mtemp3;
    if (abs(a) > abs(b)) {
        mtemp3 = b / a;
        dtemp1 = 1.0;
        dtemp2 = mtemp3.get_d();
    }
    if (abs(a) < abs(b)) {
        mtemp3 = a / b;
        dtemp1 = mtemp3.get_d();
        dtemp2 = 1.0;
    }
    if (abs(a) == abs(b)) {
        dtemp1 = 1.0;
        dtemp2 = 1.0;
    }
    mtemp3 = mpf_class(atan2(dtemp1, dtemp2));
    return mtemp3;
}

inline mpc_class sin(mpc_class a) {
    double dtemp1, dtemp2;
    mpf_class mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    dtemp1 = mtemp1.get_d();
    dtemp2 = mtemp2.get_d();
    mtemp1 = sin(dtemp1) * cosh(dtemp2);
    mtemp2 = cos(dtemp1) * sinh(dtemp2);
    mpc_class b = mpc_class(mtemp1, mtemp2);
    return b;
}

inline mpf_class log2(mpf_class x) {
    double d;
    double ln2_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln2_app = (double)exp + log10(d) / log10(2);
    return ln2_app;
}

inline mpf_class log(mpf_class x) {
    double d;
    double ln_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln_app = (double)exp * log(2.0) + log(d);
    return ln_app;
}

inline mpf_class log10(mpf_class x) {
    double d;
    double ln10_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln10_app = (double)exp * log10(2.0) + log10(d);
    return ln10_app;
}

inline mpf_class pow(mpf_class x, mplapackint y) {
    mpf_class mtemp1, mtemp2;
    if (y >= 0) {
        mpf_pow_ui(mtemp1.get_mpf_t(), x.get_mpf_t(), y);
    } else {
        mpf_pow_ui(mtemp2.get_mpf_t(), x.get_mpf_t(), -y);
        mtemp1 = 1.0 / mtemp2;
    }
    return mtemp1;
}

inline mpf_class pow(mpf_class x, mpf_class y) {
    mpf_class mtemp1, mtemp2;
    mtemp1 = y * log(x);
    mtemp2 = exp(mtemp1.get_d());
    return mtemp2;
}

inline mpf_class cos(mpf_class x) {
    mpf_class mtemp1;
    mtemp1 = cos(x.get_d());
    return mtemp1;
}

inline mpf_class sin(mpf_class x) {
    mpf_class mtemp1;
    mtemp1 = sin(x.get_d());
    return mtemp1;
}

inline mpf_class exp(mpf_class x) {
    mpf_class mtemp1;
    mtemp1 = exp(x.get_d());
    return mtemp1;
}

inline mpf_class pi(mpf_class dummy) {
    mpf_class mtemp1;
    mtemp1 = M_PI; // returns pi in double(!)
    return mtemp1;
}

inline mpc_class exp(mpc_class x) {
    mpf_class ex;
    mpf_class c;
    mpf_class s;
    mpc_class ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

static inline mpf_class cabs1(const mpc_class &z) { return abs(z.real()) + abs(z.imag()); }

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

#ifndef MPLAPACK_MINMAX_MPF_CLASS_DEFINED
#define MPLAPACK_MINMAX_MPF_CLASS_DEFINED

#include <type_traits>

// Non-template overloads for mpf_class to beat std::min/std::max templates
// when both arguments are exactly mpf_class.
inline mpf_class min(const mpf_class &a, const mpf_class &b) { return (a > b) ? b : a; }
inline mpf_class max(const mpf_class &a, const mpf_class &b) { return (a < b) ? b : a; }

inline mpf_class min(const mpf_class &a, const mpf_class &b, const mpf_class &c) {
    mpf_class r = ::min(a, b);
    return ::min(r, c);
}
inline mpf_class max(const mpf_class &a, const mpf_class &b, const mpf_class &c) {
    mpf_class r = ::max(a, b);
    return ::max(r, c);
}

// -------------------------
// 2-arg: GMP expressions
//   - Same E: this overload should win over std::max/std::min.
//   - Different E: enabled only when E1 != E2.
// -------------------------

template <class E> inline mpf_class min(const __gmp_expr<mpf_t, E> &a, const __gmp_expr<mpf_t, E> &b) {
    mpf_class aa(a), bb(b);
    return (aa > bb) ? bb : aa;
}

template <class E> inline mpf_class max(const __gmp_expr<mpf_t, E> &a, const __gmp_expr<mpf_t, E> &b) {
    mpf_class aa(a), bb(b);
    return (aa < bb) ? bb : aa;
}

template <class E1, class E2, typename = std::enable_if_t<!std::is_same_v<E1, E2>>> inline mpf_class min(const __gmp_expr<mpf_t, E1> &a, const __gmp_expr<mpf_t, E2> &b) {
    mpf_class aa(a), bb(b);
    return (aa > bb) ? bb : aa;
}

template <class E1, class E2, typename = std::enable_if_t<!std::is_same_v<E1, E2>>> inline mpf_class max(const __gmp_expr<mpf_t, E1> &a, const __gmp_expr<mpf_t, E2> &b) {
    mpf_class aa(a), bb(b);
    return (aa < bb) ? bb : aa;
}

// -------------------------
// 3-arg: Fortran semantics (NOT comparator)
//   - Third argument is any type constructible to mpf_class,
//     to prevent std::max(a,b,comp) hijack.
// -------------------------

template <class E, class C, typename = std::enable_if_t<std::is_constructible_v<mpf_class, C>>> inline mpf_class min(const __gmp_expr<mpf_t, E> &a, const __gmp_expr<mpf_t, E> &b, const C &c) {
    mpf_class r = ::min(a, b);
    return ::min(r, mpf_class(c));
}

template <class E, class C, typename = std::enable_if_t<std::is_constructible_v<mpf_class, C>>> inline mpf_class max(const __gmp_expr<mpf_t, E> &a, const __gmp_expr<mpf_t, E> &b, const C &c) {
    mpf_class r = ::max(a, b);
    return ::max(r, mpf_class(c));
}

template <class E1, class E2, class C, typename = std::enable_if_t<!std::is_same_v<E1, E2> && std::is_constructible_v<mpf_class, C>>> inline mpf_class min(const __gmp_expr<mpf_t, E1> &a, const __gmp_expr<mpf_t, E2> &b, const C &c) {
    mpf_class r = ::min(a, b);
    return ::min(r, mpf_class(c));
}

template <class E1, class E2, class C, typename = std::enable_if_t<!std::is_same_v<E1, E2> && std::is_constructible_v<mpf_class, C>>> inline mpf_class max(const __gmp_expr<mpf_t, E1> &a, const __gmp_expr<mpf_t, E2> &b, const C &c) {
    mpf_class r = ::max(a, b);
    return ::max(r, mpf_class(c));
}

// -------------------------
// 4+ args: fold (Rest must be constructible to mpf_class)
// -------------------------

template <class E, class C, class... Rest, typename = std::enable_if_t<std::is_constructible_v<mpf_class, C> && (std::is_constructible_v<mpf_class, Rest> && ...)>> inline mpf_class min(const __gmp_expr<mpf_t, E> &a, const __gmp_expr<mpf_t, E> &b, const C &c, const Rest &...rest) {
    mpf_class r = ::min(a, b, c);
    ((r = ::min(r, mpf_class(rest))), ...);
    return r;
}

template <class E, class C, class... Rest, typename = std::enable_if_t<std::is_constructible_v<mpf_class, C> && (std::is_constructible_v<mpf_class, Rest> && ...)>> inline mpf_class max(const __gmp_expr<mpf_t, E> &a, const __gmp_expr<mpf_t, E> &b, const C &c, const Rest &...rest) {
    mpf_class r = ::max(a, b, c);
    ((r = ::max(r, mpf_class(rest))), ...);
    return r;
}

template <class E1, class E2, class C, class... Rest, typename = std::enable_if_t<!std::is_same_v<E1, E2> && std::is_constructible_v<mpf_class, C> && (std::is_constructible_v<mpf_class, Rest> && ...)>> inline mpf_class min(const __gmp_expr<mpf_t, E1> &a, const __gmp_expr<mpf_t, E2> &b, const C &c, const Rest &...rest) {
    mpf_class r = ::min(a, b, c);
    ((r = ::min(r, mpf_class(rest))), ...);
    return r;
}

template <class E1, class E2, class C, class... Rest, typename = std::enable_if_t<!std::is_same_v<E1, E2> && std::is_constructible_v<mpf_class, C> && (std::is_constructible_v<mpf_class, Rest> && ...)>> inline mpf_class max(const __gmp_expr<mpf_t, E1> &a, const __gmp_expr<mpf_t, E2> &b, const C &c, const Rest &...rest) {
    mpf_class r = ::max(a, b, c);
    ((r = ::max(r, mpf_class(rest))), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_MPF_CLASS_DEFINED

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

// Integer ceil for GMP mpf_class.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_MPF_CLASS_DEFINED
#define MPLAPACK_ICEIL_MPF_CLASS_DEFINED
inline mplapackint iceil(const mpf_class &x) {
    // mpf_class -> long is trunc toward zero (via mpf_get_si).
    mplapackint t = static_cast<mplapackint>(x.get_si());
    mpf_class tt = t;
    if (x > tt) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_MPF_CLASS_DEFINED

#endif

/*
 * Copyright (c) 2008-2021
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

#ifndef _MUTILS_DD_H_
#define _MUTILS_DD_H_

#if defined ___MPLAPACK_INTERNAL___
#define DD_PRECISION 32
#define DD_PRECISION_SHORT 16

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

#include <cstring>

inline void printnum(dd_real rtmp) {
    std::cout.precision(DD_PRECISION);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

inline void printnum_short(dd_real rtmp) {
    std::cout.precision(DD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}
inline void printnum(dd_complex rtmp) {
    std::cout.precision(DD_PRECISION);
    if (rtmp.real() >= 0.0) {
        std::cout << "+" << rtmp.real();
    } else {
        std::cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
        std::cout << "+" << rtmp.imag() << "i";
    } else {
        std::cout << rtmp.imag() << "i";
    }
    return;
}
inline void printnum_short(dd_complex rtmp) {
    std::cout.precision(DD_PRECISION);
    if (rtmp.real() >= 0.0) {
        std::cout << "+" << rtmp.real();
    } else {
        std::cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
        std::cout << "+" << rtmp.imag() << "i";
    } else {
        std::cout << rtmp.imag() << "i";
    }
    return;
}
inline void sprintnum(char *buf, dd_real rtmp) {
    rtmp.write(buf, __MPLAPACK_BUFLEN__, DD_PRECISION);
    return;
}
inline void sprintnum_short(char *buf, dd_real rtmp) {
    rtmp.write(buf, __MPLAPACK_BUFLEN__, DD_PRECISION_SHORT);
    return;
}
inline void sprintnum(char *buf, dd_complex rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    rtmp.real().write(buf1, __MPLAPACK_BUFLEN__, DD_PRECISION);
    rtmp.real().write(buf2, __MPLAPACK_BUFLEN__, DD_PRECISION);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
inline void sprintnum_short(char *buf, dd_complex rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    rtmp.real().write(buf1, __MPLAPACK_BUFLEN__, DD_PRECISION_SHORT);
    rtmp.real().write(buf2, __MPLAPACK_BUFLEN__, DD_PRECISION_SHORT);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
#endif

inline dd_real pow2(dd_real a) {
    dd_real mtmp = a * a;
    return mtmp;
}

inline dd_complex pow2(dd_complex a) {
    dd_complex mtmp = a * a;
    return mtmp;
}

#include <type_traits>

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

// implementation of sign transfer function.
inline dd_real sign(dd_real a, dd_real b) {
    dd_real mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline dd_real castREAL_dd(mplapackint n) {
    dd_real ret;
    ret.x[0] = (static_cast<double>(n));
    ret.x[1] = 0.0;
    return ret;
}
inline mplapackint castINTEGER_dd(dd_real a) {
    mplapackint i = a.x[0];
    return i;
}

inline long __dd_nint(dd_real a) {
    long i;
    dd_real tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (int)tmp.x[0];
    return i;
}

inline double cast2double(dd_real a) { return a.x[0]; }

inline dd_complex sin(dd_complex a) {
    dd_real mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    dd_complex b = dd_complex(sin(mtemp1) * cosh(mtemp2), cos(mtemp1) * sinh(mtemp2));
    return b;
}

inline dd_complex cos(dd_complex a) {
    dd_real mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    dd_complex b = dd_complex(cos(mtemp1) * cosh(mtemp2), -sin(mtemp1) * sinh(mtemp2));
    return b;
}

inline dd_real log2(dd_real x) { return log10(x) / (dd_real::_log2 / dd_real::_log10); }

inline dd_complex exp(dd_complex x) {
    dd_real ex;
    dd_real c;
    dd_real s;
    dd_complex ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

inline dd_real pi(dd_real dummy) { return dd_real::_pi; }

static inline dd_real cabs1(const dd_complex &z) { return abs(z.real()) + abs(z.imag()); }

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

#ifndef MPLAPACK_MINMAX_DD_REAL_DEFINED
#define MPLAPACK_MINMAX_DD_REAL_DEFINED

inline dd_real min(const dd_real &a, const dd_real &b) { return (a > b) ? b : a; }
inline dd_real max(const dd_real &a, const dd_real &b) { return (a < b) ? b : a; }

inline dd_real min(const dd_real &a, const dd_real &b, const dd_real &c) {
    dd_real r = min(a, b);
    return min(r, c);
}
inline dd_real max(const dd_real &a, const dd_real &b, const dd_real &c) {
    dd_real r = max(a, b);
    return max(r, c);
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<dd_real, std::decay_t<Args>> && ...)>> inline dd_real min(const dd_real &a, const dd_real &b, const dd_real &c, const Args &...rest) {
    dd_real r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<dd_real, std::decay_t<Args>> && ...)>> inline dd_real max(const dd_real &a, const dd_real &b, const dd_real &c, const Args &...rest) {
    dd_real r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_DD_REAL_DEFINED

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

// Integer ceil for dd_real.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_DD_REAL_DEFINED
#define MPLAPACK_ICEIL_DD_REAL_DEFINED
inline mplapackint iceil(const dd_real &x) {
    mplapackint t = static_cast<mplapackint>(x); // trunc toward zero (via conversion)
  if (x > dd_real(t)) {
        ++t;
  }
    return t;
}
#endif // MPLAPACK_ICEIL_DD_REAL_DEFINED

#endif

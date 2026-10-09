/*
 * Copyright (c) 2008-2026
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

#ifndef _MUTILS_TD_H_
#define _MUTILS_TD_H_

#if defined MPLAPACK_INTERNAL
#define TD_PRECISION 46
#define TD_PRECISION_SHORT 16

#if !defined MPLAPACK_BUFLEN
#define MPLAPACK_BUFLEN 1024
#endif

#include <cstring>

inline void printnum(td_real rtmp) {
    std::cout.precision(TD_PRECISION);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

inline void printnum_short(td_real rtmp) {
    std::cout.precision(TD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}
inline void printnum(td_complex rtmp) {
    std::cout.precision(TD_PRECISION);
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
inline void printnum_short(td_complex rtmp) {
    std::cout.precision(TD_PRECISION);
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
inline void sprintnum(char *buf, td_real rtmp) {
    rtmp.write(buf, MPLAPACK_BUFLEN, TD_PRECISION);
    return;
}
inline void sprintnum_short(char *buf, td_real rtmp) {
    rtmp.write(buf, MPLAPACK_BUFLEN, TD_PRECISION_SHORT);
    return;
}
inline void sprintnum(char *buf, td_complex rtmp) {
    char buf1[MPLAPACK_BUFLEN], buf2[MPLAPACK_BUFLEN];
    rtmp.real().write(buf1, MPLAPACK_BUFLEN, TD_PRECISION);
    rtmp.imag().write(buf2, MPLAPACK_BUFLEN, TD_PRECISION);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
inline void sprintnum_short(char *buf, td_complex rtmp) {
    char buf1[MPLAPACK_BUFLEN], buf2[MPLAPACK_BUFLEN];
    rtmp.real().write(buf1, MPLAPACK_BUFLEN, TD_PRECISION_SHORT);
    rtmp.imag().write(buf2, MPLAPACK_BUFLEN, TD_PRECISION_SHORT);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}

#include <mplapack_hex_helpers.h>

inline void sprinthex_td(char *buf, size_t n, const td_real &x) {
    // Enough room for "[a b c]" with each limb formatted by
    // format_hex_double_fixedexp().
    char x0_buf[128];
    char x1_buf[128];
    char x2_buf[128];

    format_hex_double_fixedexp(x0_buf, sizeof(x0_buf), x.x[0]);
    format_hex_double_fixedexp(x1_buf, sizeof(x1_buf), x.x[1]);
    format_hex_double_fixedexp(x2_buf, sizeof(x2_buf), x.x[2]);

    // Caller should provide a sufficiently large buffer.
    snprintf(buf, n, "[%s %s %s]", x0_buf, x1_buf, x2_buf);
}

#endif

inline td_real pow2(td_real a) {
    td_real mtmp = a * a;
    return mtmp;
}

inline td_complex pow2(td_complex a) {
    td_complex mtmp = a * a;
    return mtmp;
}

inline td_real pow4(td_real a) {
    td_real mtmp = a * a * a * a;
    return mtmp;
}

inline td_complex pow4(td_complex a) {
    td_complex mtmp = a * a * a * a;
    return mtmp;
}

#include <type_traits>

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

// implementation of sign transfer function.
inline td_real sign(td_real a, td_real b) {
    td_real mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline td_real castREAL_td(mplapackint n) {
    td_real ret;
    ret.x[0] = (static_cast<double>(n));
    ret.x[1] = 0.0;
    ret.x[2] = 0.0;
    return ret;
}
inline mplapackint castINTEGER_td(td_real a) {
    mplapackint i = a.x[0];
    return i;
}

inline long mplapack_td_nint(td_real a) {
    long i;
    td_real tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (int)tmp.x[0];
    return i;
}

inline double cast2double(td_real a) { return a.x[0]; }

inline td_complex sin(td_complex a) {
    td_real mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    td_complex b = td_complex(sin(mtemp1) * cosh(mtemp2), cos(mtemp1) * sinh(mtemp2));
    return b;
}

inline td_complex cos(td_complex a) {
    td_real mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    td_complex b = td_complex(cos(mtemp1) * cosh(mtemp2), -sin(mtemp1) * sinh(mtemp2));
    return b;
}


inline td_complex exp(td_complex x) {
    td_real ex;
    td_real c;
    td_real s;
    td_complex ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

inline td_real pi(td_real dummy) { return td_real::_pi; }

static inline td_real cabs1(const td_complex &z) { return abs(z.real()) + abs(z.imag()); }

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

#ifndef MPLAPACK_MINMAX_TD_REAL_DEFINED
#define MPLAPACK_MINMAX_TD_REAL_DEFINED

inline td_real min(const td_real &a, const td_real &b) { return (a > b) ? b : a; }
inline td_real max(const td_real &a, const td_real &b) { return (a < b) ? b : a; }

inline td_real min(const td_real &a, const td_real &b, const td_real &c) {
    td_real r = min(a, b);
    return min(r, c);
}
inline td_real max(const td_real &a, const td_real &b, const td_real &c) {
    td_real r = max(a, b);
    return max(r, c);
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<td_real, std::decay_t<Args>> && ...)>> inline td_real min(const td_real &a, const td_real &b, const td_real &c, const Args &...rest) {
    td_real r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<td_real, std::decay_t<Args>> && ...)>> inline td_real max(const td_real &a, const td_real &b, const td_real &c, const Args &...rest) {
    td_real r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_TD_REAL_DEFINED

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

// Integer ceil for td_real.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_TD_REAL_DEFINED
#define MPLAPACK_ICEIL_TD_REAL_DEFINED
inline mplapackint iceil(const td_real &x) {
    // Truncate toward zero using the leading component.
    mplapackint t = static_cast<mplapackint>(x.x[0]);

    // Avoid ambiguous overload between td_real(int) and td_real(double).
    if (x > td_real(static_cast<double>(t))) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_TD_REAL_DEFINED

#ifndef MPLAPACK_MOD_UTILS_H
#define MPLAPACK_MOD_UTILS_H
inline mplapackint mod(mplapackint a, mplapackint b) { return a % b; }
#endif // MPLAPACK_MOD_UTILS_H

#ifndef nint
#define nint mplapack_td_nint
#endif

#endif

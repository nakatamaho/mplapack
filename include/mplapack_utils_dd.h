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

#ifndef _MUTILS_DD_H_
#define _MUTILS_DD_H_

#if defined MPLAPACK_INTERNAL
#define DD_PRECISION 30
#define DD_PRECISION_SHORT 16

#if !defined MPLAPACK_BUFLEN
#define MPLAPACK_BUFLEN 1024
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
    rtmp.write(buf, MPLAPACK_BUFLEN, DD_PRECISION);
    return;
}
inline void sprintnum_short(char *buf, dd_real rtmp) {
    rtmp.write(buf, MPLAPACK_BUFLEN, DD_PRECISION_SHORT);
    return;
}
inline void sprintnum(char *buf, dd_complex rtmp) {
    char buf1[MPLAPACK_BUFLEN], buf2[MPLAPACK_BUFLEN];
    rtmp.real().write(buf1, MPLAPACK_BUFLEN, DD_PRECISION);
    rtmp.imag().write(buf2, MPLAPACK_BUFLEN, DD_PRECISION);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
inline void sprintnum_short(char *buf, dd_complex rtmp) {
    char buf1[MPLAPACK_BUFLEN], buf2[MPLAPACK_BUFLEN];
    rtmp.real().write(buf1, MPLAPACK_BUFLEN, DD_PRECISION_SHORT);
    rtmp.imag().write(buf2, MPLAPACK_BUFLEN, DD_PRECISION_SHORT);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}

#include <mplapack_hex_helpers.h>

inline void sprinthex_dd(char *buf, size_t n, const dd_real &x) {
    // 128 is safe here because our format_hex_double_fixedexp now uses a 64-byte source
    char hi_buf[128];
    char lo_buf[128];

    format_hex_double_fixedexp(hi_buf, sizeof(hi_buf), x.x[0]);
    format_hex_double_fixedexp(lo_buf, sizeof(lo_buf), x.x[1]);

    // Ensure the output buffer 'buf' is large enough for "[hi lo]"
    // If n is at least 260, this will never truncate.
    snprintf(buf, n, "[%s %s]", hi_buf, lo_buf);
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

inline dd_real pow4(dd_real a) {
    dd_real mtmp = a * a * a * a;
    return mtmp;
}

inline dd_complex pow4(dd_complex a) {
    dd_complex mtmp = a * a * a * a;
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

inline long mplapack_dd_nint(dd_real a) {
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

// Integer ceil for dd_real.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_DD_REAL_DEFINED
#define MPLAPACK_ICEIL_DD_REAL_DEFINED
inline mplapackint iceil(const dd_real &x) {
    // Truncate toward zero using the leading component.
    mplapackint t = static_cast<mplapackint>(x.x[0]);

    // Avoid ambiguous overload between dd_real(int) and dd_real(double).
    if (x > dd_real(static_cast<double>(t))) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_DD_REAL_DEFINED

#ifndef MPLAPACK_MOD_UTILS_H
#define MPLAPACK_MOD_UTILS_H
inline mplapackint mod(mplapackint a, mplapackint b) { return a % b; }
#endif // MPLAPACK_MOD_UTILS_H

#ifndef nint
#define nint mplapack_dd_nint
#endif

#endif

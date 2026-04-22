/*
 * Copyright (c) 2008-2026
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

#include "mplapack_gmp_transcendents.h"
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

inline void sprinthex_gmp_fixed_raw(char *buf, size_t n, mpf_srcptr x) {
    if (n == 0) {
        return;
    }

    if (mpf_sgn(x) == 0) {
        const mp_bitcnt_t prec = mpf_get_prec(x);
        const size_t frac_hex_digits = static_cast<size_t>((prec - 1 + 3) / 4);
        std::string frac(frac_hex_digits, '0');
        std::snprintf(buf, n, "+0x0.%sp+0000", frac.c_str());
        return;
    }

    const bool negative = (mpf_sgn(x) < 0);
    const mp_bitcnt_t prec = mpf_get_prec(x);
    const size_t frac_hex_digits = static_cast<size_t>((prec - 1 + 3) / 4);

    mp_exp_t exp16 = 0;
    char *mant_raw = mpf_get_str(nullptr, &exp16, 16, 0, x);
    if (mant_raw == nullptr) {
        std::snprintf(buf, n, "<mpf_get_str failed>");
        return;
    }

    const char *digits = mant_raw;
    if (digits[0] == '-') {
        digits++;
    }

    if (digits[0] == '\0') {
        void (*freefunc)(void *, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefunc);
        freefunc(mant_raw, std::strlen(mant_raw) + 1);

        std::string frac(frac_hex_digits, '0');
        std::snprintf(buf, n, "%s0x0.%sp+0000", negative ? "-" : "+", frac.c_str());
        return;
    }

    auto hex_value = [](char c) -> int {
        if (c >= '0' && c <= '9') {
            return c - '0';
        }
        if (c >= 'a' && c <= 'f') {
            return 10 + (c - 'a');
        }
        if (c >= 'A' && c <= 'F') {
            return 10 + (c - 'A');
        }
        return 0;
    };

    auto hex_char = [](int v) -> char {
        static const char table[] = "0123456789abcdef";
        return table[v & 0xf];
    };

    std::string bitstream;
    bitstream.reserve(std::strlen(digits) * 4);
    for (size_t i = 0; digits[i] != '\0'; ++i) {
        const int v = hex_value(digits[i]);
        bitstream.push_back((v & 8) ? '1' : '0');
        bitstream.push_back((v & 4) ? '1' : '0');
        bitstream.push_back((v & 2) ? '1' : '0');
        bitstream.push_back((v & 1) ? '1' : '0');
    }

    const size_t first_one_pos = bitstream.find('1');
    if (first_one_pos == std::string::npos) {
        void (*freefunc)(void *, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefunc);
        freefunc(mant_raw, std::strlen(mant_raw) + 1);

        std::string frac(frac_hex_digits, '0');
        std::snprintf(buf, n, "%s0x0.%sp+0000", negative ? "-" : "+", frac.c_str());
        return;
    }

    // mpf_get_str gives value = 0.<digits> * 16^exp16
    // Binary exponent of the leading 1 is therefore:
    //   exp2 = 4*exp16 - (first_one_pos + 1)
    const long exp2 = static_cast<long>(4 * exp16 - static_cast<mp_exp_t>(first_one_pos + 1));

    std::string frac_bits = bitstream.substr(first_one_pos + 1);
    const size_t needed_bits = frac_hex_digits * 4;
    if (frac_bits.size() < needed_bits) {
        frac_bits.append(needed_bits - frac_bits.size(), '0');
    } else if (frac_bits.size() > needed_bits) {
        frac_bits.resize(needed_bits);
    }

    std::string frac;
    frac.reserve(frac_hex_digits);
    for (size_t i = 0; i < needed_bits; i += 4) {
        int v = 0;
        if (frac_bits[i + 0] == '1')
            v |= 8;
        if (frac_bits[i + 1] == '1')
            v |= 4;
        if (frac_bits[i + 2] == '1')
            v |= 2;
        if (frac_bits[i + 3] == '1')
            v |= 1;
        frac.push_back(hex_char(v));
    }

    std::snprintf(buf, n, "%s0x1.%sp%+05ld", negative ? "-" : "+", frac.c_str(), exp2);

    void (*freefunc)(void *, size_t);
    mp_get_memory_functions(nullptr, nullptr, &freefunc);
    freefunc(mant_raw, std::strlen(mant_raw) + 1);
}

inline void sprinthex_gmp_fixed(char *buf, size_t n, const mpf_class &x) { sprinthex_gmp_fixed_raw(buf, n, x.get_mpf_t()); }

#endif

inline mpf_class abs(mpf_class a) {
    mpf_class r;
    mpf_abs(r.get_mpf_t(), a.get_mpf_t());
    return r;
}

inline mpf_class sqrt(mpf_class a) {
    mpf_class r;
    mpf_sqrt(r.get_mpf_t(), a.get_mpf_t());
    return r;
}

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

inline mpf_class atan2(mpf_class a, mpf_class b) {
    return mplapack_gmp_transcendents::compute_atan2(a, b, std::max(a.get_prec(), b.get_prec()));
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
    const mp_bitcnt_t precision = x.get_prec();
    return mplapack_gmp_transcendents::div(mplapack_gmp_transcendents::compute_log(x, precision), mplapack_gmp_transcendents::log_two(precision), precision);
}

inline mpf_class log(mpf_class x) {
    return mplapack_gmp_transcendents::compute_log(x, x.get_prec());
}

inline mpf_class log10(mpf_class x) {
    const mp_bitcnt_t precision = x.get_prec();
    return mplapack_gmp_transcendents::div(mplapack_gmp_transcendents::compute_log(x, precision), mplapack_gmp_transcendents::compute_log(mpf_class(10, precision), precision), precision);
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
    return mplapack_gmp_transcendents::compute_pow(x, y, std::max(x.get_prec(), y.get_prec()));
}

inline mpf_class cos(mpf_class x) {
    return mplapack_gmp_transcendents::compute_cos(x, x.get_prec());
}

inline mpf_class sin(mpf_class x) {
    return mplapack_gmp_transcendents::compute_sin(x, x.get_prec());
}

inline mpf_class exp(mpf_class x) {
    return mplapack_gmp_transcendents::compute_exp(x, x.get_prec());
}

inline mpf_class pi(mpf_class dummy) {
    return mplapack_gmp_transcendents::pi(dummy.get_prec());
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

#ifndef MPLAPACK_MOD_UTILS_H
#define MPLAPACK_MOD_UTILS_H
inline mplapackint mod(mplapackint a, mplapackint b) { return a % b; }
#endif // MPLAPACK_MOD_UTILS_H

#endif

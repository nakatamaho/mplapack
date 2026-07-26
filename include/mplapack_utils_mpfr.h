/*
 * Copyright (c) 2008-2026
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

#include <mpfrxx_mkII.h>
#include <mpcxx_mkII.h>
using namespace mpfrxx;


#if defined ___MPLAPACK_INTERNAL___
#include <mplapack_print_double.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(mpfr_class rtmp) {
    mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp));
    return;
}

inline void printnum_short(mpfr_class rtmp) {
    mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void printnum(mpc_class ctmp) {
    mpfr_class cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_FORMAT MPFR_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

inline void printnum_short(mpc_class ctmp) {
    mpfr_class cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

inline void sprintnum(char *buf, mpfr_class rtmp) {
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum_short(char *buf, mpfr_class rtmp) {
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_SHORT_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum(char *buf, mpc_class ctmp) {
    mpfr_class cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_FORMAT MPFR_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
inline void sprintnum_short(char *buf, mpc_class ctmp) {
    mpfr_class cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

inline void sprinthex_mpfr_fixed_raw(char *buf, size_t n, mpfr_ptr x) {
    if (n == 0) {
        return;
    }

    if (mpfr_nan_p(x)) {
        std::snprintf(buf, n, "@NaN@");
        return;
    }
    if (mpfr_inf_p(x)) {
        std::snprintf(buf, n, "%s@Inf@", mpfr_sgn(x) < 0 ? "-" : "+");
        return;
    }

    const mpfr_prec_t prec = mpfr_get_prec(x);
    const size_t frac_hex_digits = static_cast<size_t>((prec - 1 + 3) / 4);

    if (mpfr_zero_p(x)) {
        std::string frac(frac_hex_digits, '0');
        std::snprintf(buf, n, "%s0x0.%sp+0000", mpfr_sgn(x) < 0 ? "-" : "+", frac.c_str());
        return;
    }

    mpfr_exp_t exp16 = 0;
    char *mant_raw = mpfr_get_str(nullptr, &exp16, 16, 0, x, MPFR_RNDN);
    if (mant_raw == nullptr) {
        std::snprintf(buf, n, "<mpfr_get_str failed>");
        return;
    }

    const bool negative = (mant_raw[0] == '-');
    const char *digits = negative ? mant_raw + 1 : mant_raw;
    if (digits[0] == '\0') {
        mpfr_free_str(mant_raw);
        std::snprintf(buf, n, "<empty mantissa>");
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

    std::string all_digits(digits);
    const int first_nibble = hex_value(all_digits[0]);

    int lead_shift = 0;
    if (first_nibble >= 8) {
        lead_shift = 3;
    } else if (first_nibble >= 4) {
        lead_shift = 2;
    } else if (first_nibble >= 2) {
        lead_shift = 1;
    } else {
        lead_shift = 0;
    }

    const long exp2 = static_cast<long>(4 * (exp16 - 1) + lead_shift);

    std::string bitstream;
    bitstream.reserve(all_digits.size() * 4 + 4);
    for (size_t i = 0; i < all_digits.size(); ++i) {
        const int v = hex_value(all_digits[i]);
        bitstream.push_back((v & 8) ? '1' : '0');
        bitstream.push_back((v & 4) ? '1' : '0');
        bitstream.push_back((v & 2) ? '1' : '0');
        bitstream.push_back((v & 1) ? '1' : '0');
    }

    const size_t first_one_pos = bitstream.find('1');
    if (first_one_pos == std::string::npos) {
        mpfr_free_str(mant_raw);
        std::string frac(frac_hex_digits, '0');
        std::snprintf(buf, n, "%s0x0.%sp+0000", negative ? "-" : "+", frac.c_str());
        return;
    }

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

    mpfr_free_str(mant_raw);
}

inline void sprinthex_mpfr_fixed(char *buf, size_t n, const mpfrxx::mpfr_class &x) {
    mpfr_ptr px = const_cast<mpfrxx::mpfr_class &>(x);
    sprinthex_mpfr_fixed_raw(buf, n, px);
}

#endif

inline mpfr_class pow2(mpfr_class a) {
    mpfr_class mtmp = a * a;
    return mtmp;
}

inline mpc_class pow2(mpc_class a) {
    mpc_class mtmp = a * a;
    return mtmp;
}

inline mpfr_class pow4(mpfr_class a) {
    mpfr_class mtmp = a * a * a * a;
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

// implementation of sign transfer function.
inline mpfr_class sign(mpfr_class a, mpfr_class b) {
    mpfr_class mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline mplapackint nint(mpfr_class a) {
    mplapackint i;
    mpfr_class tmp;
    a = a + 0.5;
    tmp = floorl(a);
    i = tmp; // cast to long
    return i;
}

inline mplapackint castINTEGER_mpfr(mpfr_class a) {
    mplapackint i;
    i = a;
    return i;
}

inline mpfr_class castREAL_mpfr(mplapackint a) {
    mpfr_class i = a;
    return i;
}

inline mpfr_class pi(mpfr_class dummy) {
    mpfr_class _PI;
    _PI = const_pi(mpfrxx::mpfr_class::default_prec);
    return _PI;
}

static inline mpfr_class cabs1(const mpc_class &z) { return abs(z.real()) + abs(z.imag()); }

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
inline mpfr_class min(const mpfr_class &a, const mpfr_class &b, const mpfr_class &c) {
    mpfr_class r = (b < a) ? b : a;
    return (c < r) ? c : r;
}
inline mpfr_class max(const mpfr_class &a, const mpfr_class &b, const mpfr_class &c) {
    mpfr_class r = (a < b) ? b : a;
    return (r < c) ? c : r;
}

// 4+ args: fold expression, mpfr_class only.
template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mpfr_class, std::decay_t<Args>> && ...)>> inline mpfr_class min(const mpfr_class &a, const mpfr_class &b, const mpfr_class &c, const Args &...rest) {
    mpfr_class r = min(a, b, c);
    ((r = (rest < r) ? rest : r), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<mpfr_class, std::decay_t<Args>> && ...)>> inline mpfr_class max(const mpfr_class &a, const mpfr_class &b, const mpfr_class &c, const Args &...rest) {
    mpfr_class r = max(a, b, c);
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

// Integer ceil for MPFR mpfr_class.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_MPREAL_DEFINED
#define MPLAPACK_ICEIL_MPREAL_DEFINED
inline mplapackint iceil(const mpfr_class &x) {
    // mpfr_class -> integer cast truncates toward zero.
    mplapackint t = static_cast<mplapackint>(x);
    if (x > mpfr_class(t)) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_MPREAL_DEFINED

#ifndef MPLAPACK_MOD_UTILS_H
#define MPLAPACK_MOD_UTILS_H
inline mplapackint mod(mplapackint a, mplapackint b) { return a % b; }
#endif // MPLAPACK_MOD_UTILS_H

#endif

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

#include <mplapack_gmpfrxx_mkII_config.h>
#include <gmpxx_mkII.h>
using namespace gmpxx;

#if defined MPLAPACK_INTERNAL
#define GMP_FORMAT "%+68.64Fe"
#define GMP_SHORT_FORMAT "%+20.16Fe"

#if !defined MPLAPACK_BUFLEN
#define MPLAPACK_BUFLEN 1024
#endif

inline void printnum(mpf_class rtmp) {
    gmp_printf(GMP_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void printnum_short(mpf_class rtmp) {
    gmp_printf(GMP_SHORT_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void printnum(mpfc_class ctmp) {
    gmp_printf(GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void printnum_short(mpfc_class ctmp) {
    gmp_printf(GMP_SHORT_FORMAT GMP_SHORT_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void sprintnum(char *buf, mpf_class rtmp) {
    gmp_snprintf(buf, MPLAPACK_BUFLEN, GMP_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void sprintnum_short(char *buf, mpf_class rtmp) {
    gmp_snprintf(buf, MPLAPACK_BUFLEN, GMP_SHORT_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void sprintnum(char *buf, mpfc_class ctmp) {
    gmp_snprintf(buf, MPLAPACK_BUFLEN, GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void sprintnum_short(char *buf, mpfc_class ctmp) {
    gmp_snprintf(buf, MPLAPACK_BUFLEN, GMP_SHORT_FORMAT GMP_SHORT_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
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



inline mpf_class pow2(mpf_class a) {
    mpf_class mtmp = a * a;
    return mtmp;
}

inline mpfc_class pow2(mpfc_class a) {
    mpfc_class mtmp = a * a;
    return mtmp;
}

inline mpf_class pow4(mpf_class a) {
    mpf_class mtmp = a * a * a * a;
    return mtmp;
}

inline mpfc_class pow4(mpfc_class a) {
    mpfc_class mtmp = a * a * a * a;
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
    return a.get_integer<mplapackint>();
}

inline mplapackint nint(mpf_class a) {
    mpf_class tmp;
    a = a + 0.5;
    mpf_floor(tmp.get_mpf_t(), a.get_mpf_t());
    return tmp.get_integer<mplapackint>();
}

inline double cast2double(mpf_class a) { return a.get_d(); }

inline mpf_class pi(mpf_class dummy) { return gmpxx::pi(dummy.get_prec()); }
inline mpf_class e(mpf_class dummy) { return gmpxx::e(dummy.get_prec()); }
inline mpf_class log_ten(mpf_class dummy) { return gmpxx::log_ten(dummy.get_prec()); }
inline mpf_class inv_log_two(mpf_class dummy) { return gmpxx::inv_log_two(dummy.get_prec()); }
inline mpf_class pi_over_two(mpf_class dummy) { return gmpxx::pi_over_two(dummy.get_prec()); }
inline mpf_class pi_over_four(mpf_class dummy) { return gmpxx::pi_over_four(dummy.get_prec()); }
inline mpf_class two_pi(mpf_class dummy) { return gmpxx::two_pi(dummy.get_prec()); }

static inline mpf_class cabs1(const mpfc_class &z) { return gmpxx::abs(z.real()) + gmpxx::abs(z.imag()); }

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

template <typename T>
inline constexpr bool mplapack_gmp_mpf_operand_v =
    gmpfrxx_mkII::detail::is_mpf_expression_operand_v<T>;

template <typename... Ts>
inline constexpr bool mplapack_gmp_has_mpf_object_or_node_v =
    (gmpfrxx_mkII::detail::is_mpf_object_or_node_v<Ts> || ...);

template <typename First, typename... Rest>
inline mpf_class mplapack_gmp_min_values(const First &first, const Rest &...rest) {
    mpf_class result(first);
    auto update = [&](const auto &value) {
        mpf_class candidate(value);
        if (result > candidate) {
            result = candidate;
        }
    };
    (update(rest), ...);
    return result;
}

template <typename First, typename... Rest>
inline mpf_class mplapack_gmp_max_values(const First &first, const Rest &...rest) {
    mpf_class result(first);
    auto update = [&](const auto &value) {
        mpf_class candidate(value);
        if (result < candidate) {
            result = candidate;
        }
    };
    (update(rest), ...);
    return result;
}

inline mpf_class min(const mpf_class &a, const mpf_class &b) { return mplapack_gmp_min_values(a, b); }
inline mpf_class max(const mpf_class &a, const mpf_class &b) { return mplapack_gmp_max_values(a, b); }

inline mpf_class min(const mpf_class &a, const mpf_class &b, const mpf_class &c) { return mplapack_gmp_min_values(a, b, c); }
inline mpf_class max(const mpf_class &a, const mpf_class &b, const mpf_class &c) { return mplapack_gmp_max_values(a, b, c); }

template <typename A, typename B,
          std::enable_if_t<mplapack_gmp_mpf_operand_v<A> &&
                               mplapack_gmp_mpf_operand_v<B> &&
                               mplapack_gmp_has_mpf_object_or_node_v<A, B>,
                           int> = 0>
inline mpf_class min(const A &a, const B &b) { return mplapack_gmp_min_values(a, b); }

template <typename A, typename B,
          std::enable_if_t<mplapack_gmp_mpf_operand_v<A> &&
                               mplapack_gmp_mpf_operand_v<B> &&
                               mplapack_gmp_has_mpf_object_or_node_v<A, B>,
                           int> = 0>
inline mpf_class max(const A &a, const B &b) { return mplapack_gmp_max_values(a, b); }

template <typename A, typename B, typename C,
          std::enable_if_t<mplapack_gmp_mpf_operand_v<A> &&
                               mplapack_gmp_mpf_operand_v<B> &&
                               mplapack_gmp_mpf_operand_v<C> &&
                               mplapack_gmp_has_mpf_object_or_node_v<A, B, C>,
                           int> = 0>
inline mpf_class min(const A &a, const B &b, const C &c) { return mplapack_gmp_min_values(a, b, c); }

template <typename A, typename B, typename C,
          std::enable_if_t<mplapack_gmp_mpf_operand_v<A> &&
                               mplapack_gmp_mpf_operand_v<B> &&
                               mplapack_gmp_mpf_operand_v<C> &&
                               mplapack_gmp_has_mpf_object_or_node_v<A, B, C>,
                           int> = 0>
inline mpf_class max(const A &a, const B &b, const C &c) { return mplapack_gmp_max_values(a, b, c); }

template <typename A, typename B, typename C, typename D, typename... Rest,
          std::enable_if_t<mplapack_gmp_mpf_operand_v<A> &&
                               mplapack_gmp_mpf_operand_v<B> &&
                               mplapack_gmp_mpf_operand_v<C> &&
                               mplapack_gmp_mpf_operand_v<D> &&
                               (mplapack_gmp_mpf_operand_v<Rest> && ...) &&
                               mplapack_gmp_has_mpf_object_or_node_v<A, B, C, D, Rest...>,
                           int> = 0>
inline mpf_class min(const A &a, const B &b, const C &c, const D &d, const Rest &...rest) {
    return mplapack_gmp_min_values(a, b, c, d, rest...);
}

template <typename A, typename B, typename C, typename D, typename... Rest,
          std::enable_if_t<mplapack_gmp_mpf_operand_v<A> &&
                               mplapack_gmp_mpf_operand_v<B> &&
                               mplapack_gmp_mpf_operand_v<C> &&
                               mplapack_gmp_mpf_operand_v<D> &&
                               (mplapack_gmp_mpf_operand_v<Rest> && ...) &&
                               mplapack_gmp_has_mpf_object_or_node_v<A, B, C, D, Rest...>,
                           int> = 0>
inline mpf_class max(const A &a, const B &b, const C &c, const D &d, const Rest &...rest) {
    return mplapack_gmp_max_values(a, b, c, d, rest...);
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
    // mpf_class -> integer truncates toward zero.
    mplapackint t = x.get_integer<mplapackint>();
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

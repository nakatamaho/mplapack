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

#ifndef MPLAPACK_HAVE_STD_ABS_FLOAT64X
#define MPLAPACK_HAVE_STD_ABS_FLOAT64X 0
#endif
#ifndef MPLAPACK_HAVE_STD_MATH_FLOAT64X
#define MPLAPACK_HAVE_STD_MATH_FLOAT64X 0
#endif

#if defined ___MPLAPACK_INTERNAL___

#include <cstring>
inline uint16_t load_u16_be(const unsigned char *p) { return (static_cast<uint16_t>(p[0]) << 8) | (static_cast<uint16_t>(p[1]) << 0); }
inline uint16_t load_u16_le(const unsigned char *p) { return (static_cast<uint16_t>(p[1]) << 8) | (static_cast<uint16_t>(p[0]) << 0); }
inline uint64_t load_u64_be(const unsigned char *p) { return (static_cast<uint64_t>(p[0]) << 56) | (static_cast<uint64_t>(p[1]) << 48) | (static_cast<uint64_t>(p[2]) << 40) | (static_cast<uint64_t>(p[3]) << 32) | (static_cast<uint64_t>(p[4]) << 24) | (static_cast<uint64_t>(p[5]) << 16) | (static_cast<uint64_t>(p[6]) << 8) | (static_cast<uint64_t>(p[7]) << 0); }
inline uint64_t load_u64_le(const unsigned char *p) { return (static_cast<uint64_t>(p[7]) << 56) | (static_cast<uint64_t>(p[6]) << 48) | (static_cast<uint64_t>(p[5]) << 40) | (static_cast<uint64_t>(p[4]) << 32) | (static_cast<uint64_t>(p[3]) << 24) | (static_cast<uint64_t>(p[2]) << 16) | (static_cast<uint64_t>(p[1]) << 8) | (static_cast<uint64_t>(p[0]) << 0); }

inline void sprinthex_binary80_bits_le(char *buf, size_t n, const unsigned char bytes[10]) {
    if (n == 0) {
        return;
    }

    // x87 little-endian memory layout:
    // bytes[0..7]  = significand, little-endian
    // bytes[8..9]  = sign/exponent, little-endian
    const uint64_t signif = load_u64_le(bytes + 0);
    const uint16_t se = load_u16_le(bytes + 8);

    const bool negative = ((se >> 15) & 1u) != 0;
    const uint16_t exp_raw = static_cast<uint16_t>(se & 0x7fffu);

    const bool int_bit = ((signif >> 63) & 1u) != 0;
    const uint64_t frac63 = (signif & 0x7fffffffffffffffULL);

    const char *sign_str = negative ? "-" : "+";

    // Zero
    if (exp_raw == 0 && !int_bit && frac63 == 0) {
        std::snprintf(buf, n, "%s0x0.0000000000000000p+0000", sign_str);
        return;
    }

    // Inf / NaN / unsupported pseudo forms
    if (exp_raw == 0x7fff) {
        if (int_bit && frac63 == 0) {
            std::snprintf(buf, n, "%s@Inf@", sign_str);
        } else {
            std::snprintf(buf, n, "@NaN@");
        }
        return;
    }

    auto hex_digit = [](int v) -> char {
        static const char table[] = "0123456789abcdef";
        return table[v & 0xf];
    };

    // Normal numbers in x87 extended have explicit integer bit = 1.
    if (exp_raw != 0 && int_bit) {
        // Fraction is 63 bits, so fixed hex fraction uses 16 digits = 64 bits,
        // with the last low bit padded by 0.
        char frac_hex[17];
        for (int i = 0; i < 16; ++i) {
            int v = 0;
            for (int j = 0; j < 4; ++j) {
                const int bit_index = i * 4 + j; // 0..63 over displayed fraction bits
                int bit = 0;
                if (bit_index < 63) {
                    const int src_bit = 62 - bit_index;
                    bit = static_cast<int>((frac63 >> src_bit) & 1ULL);
                } else {
                    bit = 0; // pad
                }
                if (bit) {
                    v |= (1 << (3 - j));
                }
            }
            frac_hex[i] = hex_digit(v);
        }
        frac_hex[16] = '\0';

        const int exp_unbiased = static_cast<int>(exp_raw) - 16383;
        std::snprintf(buf, n, "%s0x1.%sp%+06d", sign_str, frac_hex, exp_unbiased);
        return;
    }

    // Subnormal / pseudo-denormal / unnormal:
    // normalize from raw bitstream of integer-bit + fraction bits (64 bits total).
    {
        std::string sig_bits;
        sig_bits.reserve(64);

        sig_bits.push_back(int_bit ? '1' : '0');
        for (int i = 62; i >= 0; --i) {
            sig_bits.push_back(((frac63 >> i) & 1ULL) ? '1' : '0');
        }

        const size_t first_one = sig_bits.find('1');
        if (first_one == std::string::npos) {
            std::snprintf(buf, n, "%s0x0.0000000000000000p+0000", sign_str);
            return;
        }

        // For subnormals in 80-bit extended, exponent is 1-bias for leading bit search.
        // Then shift according to first_one.
        const int exp_unbiased = -16382 - static_cast<int>(first_one);

        std::string mant_bits = sig_bits.substr(first_one + 1);
        if (mant_bits.size() < 64) {
            mant_bits.append(64 - mant_bits.size(), '0');
        } else if (mant_bits.size() > 64) {
            mant_bits.resize(64);
        }

        char frac_hex[17];
        for (int i = 0; i < 16; ++i) {
            int v = 0;
            for (int j = 0; j < 4; ++j) {
                if (mant_bits[i * 4 + j] == '1') {
                    v |= (1 << (3 - j));
                }
            }
            frac_hex[i] = hex_digit(v);
        }
        frac_hex[16] = '\0';

        std::snprintf(buf, n, "%s0x1.%sp%+06d", sign_str, frac_hex, exp_unbiased);
    }
}

inline void sprinthex_binary80_bits_be(char *buf, size_t n, const unsigned char bytes[10]) {
    if (n == 0) {
        return;
    }

    // Big-endian variant:
    // bytes[0..1] = sign/exponent, big-endian
    // bytes[2..9] = significand, big-endian
    const uint16_t se = load_u16_be(bytes + 0);
    const uint64_t signif = load_u64_be(bytes + 2);

    const bool negative = ((se >> 15) & 1u) != 0;
    const uint16_t exp_raw = static_cast<uint16_t>(se & 0x7fffu);

    const bool int_bit = ((signif >> 63) & 1u) != 0;
    const uint64_t frac63 = (signif & 0x7fffffffffffffffULL);

    const char *sign_str = negative ? "-" : "+";

    if (exp_raw == 0 && !int_bit && frac63 == 0) {
        std::snprintf(buf, n, "%s0x0.0000000000000000p+0000", sign_str);
        return;
    }

    if (exp_raw == 0x7fff) {
        if (int_bit && frac63 == 0) {
            std::snprintf(buf, n, "%s@Inf@", sign_str);
        } else {
            std::snprintf(buf, n, "@NaN@");
        }
        return;
    }

    auto hex_digit = [](int v) -> char {
        static const char table[] = "0123456789abcdef";
        return table[v & 0xf];
    };

    if (exp_raw != 0 && int_bit) {
        char frac_hex[17];
        for (int i = 0; i < 16; ++i) {
            int v = 0;
            for (int j = 0; j < 4; ++j) {
                const int bit_index = i * 4 + j;
                int bit = 0;
                if (bit_index < 63) {
                    const int src_bit = 62 - bit_index;
                    bit = static_cast<int>((frac63 >> src_bit) & 1ULL);
                } else {
                    bit = 0;
                }
                if (bit) {
                    v |= (1 << (3 - j));
                }
            }
            frac_hex[i] = hex_digit(v);
        }
        frac_hex[16] = '\0';

        const int exp_unbiased = static_cast<int>(exp_raw) - 16383;
        std::snprintf(buf, n, "%s0x1.%sp%+06d", sign_str, frac_hex, exp_unbiased);
        return;
    }

    std::string sig_bits;
    sig_bits.reserve(64);

    sig_bits.push_back(int_bit ? '1' : '0');
    for (int i = 62; i >= 0; --i) {
        sig_bits.push_back(((frac63 >> i) & 1ULL) ? '1' : '0');
    }

    const size_t first_one = sig_bits.find('1');
    if (first_one == std::string::npos) {
        std::snprintf(buf, n, "%s0x0.0000000000000000p+0000", sign_str);
        return;
    }

    const int exp_unbiased = -16382 - static_cast<int>(first_one);

    std::string mant_bits = sig_bits.substr(first_one + 1);
    if (mant_bits.size() < 64) {
        mant_bits.append(64 - mant_bits.size(), '0');
    } else if (mant_bits.size() > 64) {
        mant_bits.resize(64);
    }

    char frac_hex[17];
    for (int i = 0; i < 16; ++i) {
        int v = 0;
        for (int j = 0; j < 4; ++j) {
            if (mant_bits[i * 4 + j] == '1') {
                v |= (1 << (3 - j));
            }
        }
        frac_hex[i] = hex_digit(v);
    }
    frac_hex[16] = '\0';

    std::snprintf(buf, n, "%s0x1.%sp%+06d", sign_str, frac_hex, exp_unbiased);
}

template <class T> inline void sprinthex_binary80_fixed(char *buf, size_t n, const T &x) {
    static_assert(sizeof(T) >= 10, "binary80 object must be at least 10 bytes");
    unsigned char bytes[sizeof(T)];
    ::memcpy(bytes, &x, sizeof(T));

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    // On common x86/x86_64 ABIs, the first 10 bytes contain the x87 payload.
    sprinthex_binary80_bits_le(buf, n, bytes);
#else
    // Conservative big-endian assumption: first 10 bytes are the x87 payload.
    sprinthex_binary80_bits_be(buf, n, bytes);
#endif
}

#if MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL

#define LDBL_FORMAT "%+25.21Le"
#define LDBL_SHORT_FORMAT "%+20.16Le"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(long double rtmp) { printf(LDBL_FORMAT, rtmp); }
inline void printnum(std::complex<long double> ctmp) { printf(LDBL_FORMAT LDBL_FORMAT "i", ctmp.real(), ctmp.imag()); }
inline void printnum_short(long double rtmp) { printf(LDBL_SHORT_FORMAT, rtmp); }
inline void printnum_short(std::complex<long double> ctmp) { printf(LDBL_SHORT_FORMAT LDBL_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }
inline void sprintnum(char *buf, long double rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, LDBL_FORMAT, rtmp); }
inline void sprintnum(char *buf, std::complex<long double> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, LDBL_FORMAT LDBL_FORMAT "i", ctmp.real(), ctmp.imag()); }
inline void sprintnum_short(char *buf, long double rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, LDBL_SHORT_FORMAT, rtmp); }
inline void sprintnum_short(char *buf, std::complex<long double> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, LDBL_SHORT_FORMAT LDBL_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }

#elif MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_STRFROMF64X

#define FLOAT64X_FORMAT "%.21e"
#define FLOAT64X_SHORT_FORMAT "%.16e"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

// Must be before including <stdlib.h> on glibc to expose strfromf64x.
#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif

#include <stdlib.h> // strfromf64x
#include <stdio.h>  // fputs, putchar
#include <string.h> // memmove, strlen, memcpy
#include <complex>

// Helper: convert _Float64x to string using strfromf64x, then explicitly
// prepend '+' for non-negative values (including +0, +inf, +nan).
// strfromf64x does not portably support the '+' flag, so we do it manually.
// Returns number of chars written (excluding NUL), or negative on error.
static inline int mplapack_strfromf64x(char *buf, size_t buflen, const char *fmt, _Float64x x) {
    if (buf == nullptr || buflen == 0 || fmt == nullptr)
        return -1;
    // Use a staging buffer so we have room to insert '+' without clobbering buf.
    char tmp[__MPLAPACK_BUFLEN__];
    size_t tmplen = (buflen <= sizeof(tmp)) ? buflen : sizeof(tmp);
    int n = strfromf64x(tmp, tmplen, fmt, x);
    if (n < 0) {
        buf[0] = '\0';
        return n;
    }
    tmp[tmplen - 1] = '\0'; // always NUL-terminate
    // If the result does not start with '-', explicitly prepend '+'.
    if (tmp[0] != '-') {
        size_t len = strlen(tmp);
        if (len + 2 <= tmplen) {
            memmove(tmp + 1, tmp, len + 1); // shift right, carries NUL
            tmp[0] = '+';
            n = (int)(len + 1);
        }
        // If no room (extremely truncated), leave without '+' rather than corrupt.
    }
    // Copy into caller-provided buffer.
    size_t copy_len = strlen(tmp);
    if (copy_len >= buflen)
        copy_len = buflen - 1;
    memcpy(buf, tmp, copy_len);
    buf[copy_len] = '\0';
    return (int)copy_len;
}

// printnum / printnum_short: stringize then print
inline void printnum(_Float64x rtmp) {
    char buf[__MPLAPACK_BUFLEN__];
    mplapack_strfromf64x(buf, sizeof(buf), FLOAT64X_FORMAT, rtmp);
    fputs(buf, stdout);
}

inline void printnum(std::complex<_Float64x> ctmp) {
    char re[__MPLAPACK_BUFLEN__];
    char im[__MPLAPACK_BUFLEN__];
    mplapack_strfromf64x(re, sizeof(re), FLOAT64X_FORMAT, ctmp.real());
    mplapack_strfromf64x(im, sizeof(im), FLOAT64X_FORMAT, ctmp.imag());
    fputs(re, stdout);
    fputs(im, stdout);
    putchar('i');
}

inline void printnum_short(_Float64x rtmp) {
    char buf[__MPLAPACK_BUFLEN__];
    mplapack_strfromf64x(buf, sizeof(buf), FLOAT64X_SHORT_FORMAT, rtmp);
    fputs(buf, stdout);
}

inline void printnum_short(std::complex<_Float64x> ctmp) {
    char re[__MPLAPACK_BUFLEN__];
    char im[__MPLAPACK_BUFLEN__];
    mplapack_strfromf64x(re, sizeof(re), FLOAT64X_SHORT_FORMAT, ctmp.real());
    mplapack_strfromf64x(im, sizeof(im), FLOAT64X_SHORT_FORMAT, ctmp.imag());
    fputs(re, stdout);
    fputs(im, stdout);
    putchar('i');
}

// sprintnum / sprintnum_short: write into caller-provided buf
inline void sprintnum(char *buf, _Float64x rtmp) { mplapack_strfromf64x(buf, __MPLAPACK_BUFLEN__, FLOAT64X_FORMAT, rtmp); }
inline void sprintnum(char *buf, std::complex<_Float64x> ctmp) {
    if (buf == nullptr)
        return;
    // Write real part
    int nre = mplapack_strfromf64x(buf, __MPLAPACK_BUFLEN__, FLOAT64X_FORMAT, ctmp.real());
    if (nre < 0) {
        buf[0] = '\0';
        return;
    }
    // Append imag part into remaining space
    size_t used = 0;
    while (used + 1 < (size_t)__MPLAPACK_BUFLEN__ && buf[used] != '\0')
        used++;

    if (used + 1 >= (size_t)__MPLAPACK_BUFLEN__) {
        buf[__MPLAPACK_BUFLEN__ - 1] = '\0';
        return;
    }
    int nim = mplapack_strfromf64x(buf + used, (size_t)__MPLAPACK_BUFLEN__ - used, FLOAT64X_FORMAT, ctmp.imag());
    if (nim < 0) {
        buf[used] = '\0';
        return;
    }
    // Append 'i' if space
    used = 0;
    while (used + 1 < (size_t)__MPLAPACK_BUFLEN__ && buf[used] != '\0')
        used++;
    if (used + 2 <= (size_t)__MPLAPACK_BUFLEN__) {
        buf[used] = 'i';
        buf[used + 1] = '\0';
    } else {
        buf[__MPLAPACK_BUFLEN__ - 1] = '\0';
    }
}
inline void sprintnum_short(char *buf, _Float64x rtmp) { (void)mplapack_strfromf64x(buf, __MPLAPACK_BUFLEN__, FLOAT64X_SHORT_FORMAT, rtmp); }
inline void sprintnum_short(char *buf, std::complex<_Float64x> ctmp) {
    if (buf == nullptr)
        return;
    int nre = mplapack_strfromf64x(buf, __MPLAPACK_BUFLEN__, FLOAT64X_SHORT_FORMAT, ctmp.real());
    if (nre < 0) {
        buf[0] = '\0';
        return;
    }
    size_t used = 0;
    while (used + 1 < (size_t)__MPLAPACK_BUFLEN__ && buf[used] != '\0')
        used++;
    if (used + 1 >= (size_t)__MPLAPACK_BUFLEN__) {
        buf[__MPLAPACK_BUFLEN__ - 1] = '\0';
        return;
    }
    int nim = mplapack_strfromf64x(buf + used, (size_t)__MPLAPACK_BUFLEN__ - used, FLOAT64X_SHORT_FORMAT, ctmp.imag());
    if (nim < 0) {
        buf[used] = '\0';
        return;
    }
    used = 0;
    while (used + 1 < (size_t)__MPLAPACK_BUFLEN__ && buf[used] != '\0')
        used++;
    if (used + 2 <= (size_t)__MPLAPACK_BUFLEN__) {
        buf[used] = 'i';
        buf[used + 1] = '\0';
    } else {
        buf[__MPLAPACK_BUFLEN__ - 1] = '\0';
    }
}
#else
#error
#endif // MPLAPACK_BINARY80_IO
#endif // ___MPLAPACK_INTERNAL___

#if MPLAPACK_BINARY80_MATH == MPLAPACK_BINARY80_MATH_LDBL
#include <cmath> // do not rely on transitive includes

// Expose <cmath> overload sets for unqualified calls in auto-converted LAPACK code.
// This ensures long double stays long double (avoids accidental long double -> double demotion),
// and avoids ambiguous overloads caused by custom global wrappers.
using std::abs;
using std::atan2;
using std::cos;
using std::cosh;
using std::exp;
using std::floor;
using std::log;
using std::log10;
using std::log2;
using std::pow;
using std::sin;
using std::sinh;
using std::sqrt;

inline long double pow(const long &a, const long &b) { return powl((long double)a, (long double)b); }
inline long double pow(const int &a, const long &b) { return powl((long double)a, (long double)b); }
inline long double pow(const long double &a, const long &b) { return powl(a, (long double)b); }
inline long double pow2(const long double &a) { return a * a; }
inline long double pow4(const long double &a) { return a * a * a * a; }
inline std::complex<long double> pow2(const std::complex<long double> &a) { return a * a; }
inline std::complex<long double> pow4(const std::complex<long double> &a) { return a * a * a * a; }

inline long double ldexp(const long double &a, int exp) { return ldexpl(a, exp); }
inline long double nextafter(const long double &a, const long double &b) { return nextafterl(a, b); }

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
inline long double pi([[maybe_unused]] long double dummy) { return 0xc.90fdaa22168c235p-2L; }
#elif MPLAPACK_BINARY80_MATH == MPLAPACK_BINARY80_MATH_F64X
#pragma once
#include <cmath>
#include <math.h>
#include <complex>
// Basic real functions
#if MPLAPACK_HAVE_STD_MATH_FLOAT64X == 1
using std::atan2;
using std::cos;
using std::cosh;
using std::exp;
using std::floor;
using std::ldexp;
using std::log;
using std::log10;
using std::log2;
using std::nextafter;
using std::pow;
using std::sin;
using std::sinh;
using std::sqrt;
#else
inline _Float64x atan2(_Float64x a, _Float64x b) { return ::atan2f64x(a, b); }
inline _Float64x cos(_Float64x a) { return ::cosf64x(a); }
inline _Float64x cosh(_Float64x a) { return ::coshf64x(a); }
inline _Float64x exp(_Float64x a) { return ::expf64x(a); }
inline _Float64x floor(_Float64x a) { return ::floorf64x(a); }
inline _Float64x log(_Float64x a) { return ::logf64x(a); }
inline _Float64x log10(_Float64x a) { return ::log10f64x(a); }
inline _Float64x log2(_Float64x a) { return ::log2f64x(a); }
inline _Float64x pow(_Float64x a, _Float64x b) { return ::powf64x(a, b); }
inline _Float64x sin(_Float64x a) { return ::sinf64x(a); }
inline _Float64x sinh(_Float64x a) { return ::sinhf64x(a); }
inline _Float64x sqrt(_Float64x a) { return ::sqrtf64x(a); }
#endif

// Integer-exponent convenience overloads (matches your long double intent)
inline _Float64x pow(const long &a, const long &b) { return ::powf64x((_Float64x)a, (_Float64x)b); }
inline _Float64x pow(const int &a, const long &b) { return ::powf64x((_Float64x)a, (_Float64x)b); }
inline _Float64x pow(const _Float64x &a, const long &b) { return ::powf64x(a, (_Float64x)b); }

// pow2/pow4 (real/complex)
inline _Float64x pow2(const _Float64x &a) { return a * a; }
inline _Float64x pow4(const _Float64x &a) {
    _Float64x t = a * a;
    return t * t;
}

inline std::complex<_Float64x> pow2(const std::complex<_Float64x> &a) { return a * a; }
inline std::complex<_Float64x> pow4(const std::complex<_Float64x> &a) {
    auto t = a * a;
    return t * t;
}

// ldexp/nextafter
#if MPLAPACK_HAVE_STD_MATH_FLOAT64X != 1
inline _Float64x ldexp(const _Float64x &a, int exp) { return ::ldexpf64x(a, exp); }
inline _Float64x nextafter(const _Float64x &a, const _Float64x &b) { return ::nextafterf64x(a, b); }
#endif

// absolute value
#if MPLAPACK_HAVE_STD_ABS_FLOAT64X == 1
using std::abs;
#else
inline _Float64x abs(_Float64x a) { return ::fabsf64x(a); }
#endif

// sign transfer (your semantics preserved)
inline _Float64x sign(const _Float64x &a, const _Float64x &b) {
    _Float64x mtmp = abs(a);
    if (b < (_Float64x)0.0)
        mtmp = -mtmp;
    return mtmp;
}

// cast
inline double cast2double(_Float64x a) { return (double)a; }

// nint: your exact behavior preserved (ties go away-from-zero due to +0.5 then floor)
inline long nint(_Float64x a) {
    a = a + (_Float64x)0.5;
    _Float64x tmp = floor(a);
    long i = (long)tmp;
    return i;
}

// integer/real casts
inline _Float64x castREAL_binary80(mplapackint n) { return (_Float64x)n; }
inline mplapackint castINTEGER_binary80(_Float64x a) { return (mplapackint)a; }

// pi: avoid M_PIl dependence for _Float64x; use a hex float and cast.
inline _Float64x pi(_Float64x /*dummy*/) { return (_Float64x)0xc.90fdaa22168c235p-2L; }
#else
#error
#endif // MPLAPACK_BINARY80_MATH

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

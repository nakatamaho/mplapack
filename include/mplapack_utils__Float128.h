/*
 * Copyright (c) 2012-2025
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

#ifndef _MUTILS__FLOAT128_H_
#define _MUTILS__FLOAT128_H_

#include "mplapack_config.h"

#if (defined ___MPLAPACK__FLOAT128_ONLY___ || defined ___MPLAPACK__FLOAT128_IS_LONGDOUBLE___)

// Intel C/C++ compilers currently only supports __float128
// and following workaround works
#if defined(__INTEL_COMPILER)
#define _Float128 __float128
#endif

#if defined ___MPLAPACK_INTERNAL___
#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif
#include <string.h>

#define FLOAT128_FORMAT "%.40e"
#define FLOAT128_SHORT_FORMAT "%.16e"

inline void printnum(_Float128 rtmp) {
    char buf[__MPLAPACK_BUFLEN__];
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp);
    if (rtmp >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
}
inline void printnum_short(_Float128 rtmp) {
    char buf[__MPLAPACK_BUFLEN__];
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp);
    if (rtmp >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
}
inline void printnum(std::complex<_Float128> rtmp) {
    char buf[__MPLAPACK_BUFLEN__];
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp.real());
    if (rtmp.real() >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp.imag());
    if (rtmp.imag() >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
    printf("i");
}
inline void printnum_short(std::complex<_Float128> rtmp) {
    char buf[__MPLAPACK_BUFLEN__];
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp.real());
    if (rtmp.real() >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
    strfromf128(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp.imag());
    if (rtmp.imag() >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
    printf("i");
}
inline void sprintnum(char *buf, _Float128 rtmp) {
    char buf1[__MPLAPACK_BUFLEN__];
    buf[0] = '\0';
    if (rtmp >= 0.0)
        strcat(buf, "+");
    strfromf128(buf1, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp);
    strcat(buf, buf1);
}
inline void sprintnum_short(char *buf, _Float128 rtmp) {
    char buf1[__MPLAPACK_BUFLEN__];
    buf[0] = '\0';
    if (rtmp >= 0.0)
        strcat(buf, "+");
    strfromf128(buf1, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp);
    strcat(buf, buf1);
}
inline void sprintnum(char *buf, std::complex<_Float128> rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    buf[0] = '\0';
    if (rtmp.real() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf1, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp.real());
    strcat(buf, buf1);
    if (rtmp.imag() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf2, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp.imag());
    strcat(buf, buf2);
    strcat(buf, "i");
}
inline void sprintnum_short(char *buf, std::complex<_Float128> rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    buf[0] = '\0';
    if (rtmp.real() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf1, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp.real());
    strcat(buf, buf1);
    if (rtmp.imag() >= 0.0)
        strcat(buf, "+");
    strfromf128(buf2, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp.imag());
    strcat(buf, buf2);
    strcat(buf, "i");
}
#endif

// when _Float128 == long double, followings are already defined.

#if !defined ___MPLAPACK__FLOAT128_IS_LONGDOUBLE___

inline _Float128 pow(const _Float128 &a, const _Float128 &b) { return powf128(a, b); }
inline _Float128 pow(const long &a, const long &b) { return powf128((_Float128)a, (_Float128)b); }
inline _Float128 pow(const int &a, const long &b) { return powf128((_Float128)a, (_Float128)b); }
inline _Float128 pow(const _Float128 &a, const long &b) { return powf128(a, (_Float128)b); }
inline _Float128 sqrt(const _Float128 &a) { return sqrtf128(a); }

inline _Float128 sin(_Float128 a) { return sinf128(a); }
inline _Float128 sinh(_Float128 a) { return sinhf128(a); }
inline _Float128 cos(_Float128 a) { return cosf128(a); }
inline _Float128 cosh(_Float128 a) { return coshf128(a); }

inline _Float128 atan2(_Float128 a, _Float128 b) { return atan2f128(a, b); }

inline _Float128 exp(const _Float128 &a) { return expf128(a); }
inline _Float128 log(const _Float128 &a) { return logf128(a); }
inline _Float128 log10(const _Float128 &a) { return log10f128(a); }
inline _Float128 log2(const _Float128 &a) { return logf128(a) / logf128(2.0); }
inline _Float128 ceil(_Float128 a) { return ceilf128(a); }
inline _Float128 ldexp(const _Float128 &a, int exp) { return ldexpf128(a, exp); }
inline _Float128 nextafter(const _Float128 &a, const _Float128 &b) { return nextafterf128(a, b); }

#else
// _Float128 is long double
inline _Float128 ldexp(const _Float128 &a, int exp) { return ldexpl(a, exp); }
inline _Float128 nextafter(const _Float128 &a, const _Float128 &b) { return nextafterl(a, b); }

#endif

#ifdef __cplusplus
extern "C" {
#endif
#include <complex.h>
#include <complex>
#include <math.h>
#ifdef __cplusplus
}
#endif

inline _Float128 abs(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    _Float128 c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    c = cabsf128(b);
    return c;
}

inline std::complex<_Float128> sqrt(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csqrtf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> sin(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csinf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> cos(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = ccosf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> exp(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = cexpf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> log(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = clogf128(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

// implementation of sign transfer function.
inline _Float128 sign(_Float128 a, _Float128 b) {
    _Float128 mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline _Float128 castREAL__Float128(mplapackint n) {
    _Float128 ret = n;
    return ret;
}

inline mplapackint castINTEGER__Float128(_Float128 a) {
    mplapackint i = a;
    return i;
}

inline long nint(_Float128 a) {
    long i;
    _Float128 tmp;
    a = a + 0.5;
    tmp = floorf128(a);
    i = (long)tmp;
    return i;
}
inline double cast2double(_Float128 a) { return (double)a; }

inline _Float128 pi(_Float128 dummy) { return M_PIf128; }
#endif // (defined ___MPLAPACK__FLOAT128_ONLY___ || defined ___MPLAPACK__FLOAT128_IS_LONGDOUBLE___)

#if defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___

#if defined ___MPLAPACK_INTERNAL___
#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif
#include <string.h>

#define FLOAT128_FORMAT "%+.40Le"
#define FLOAT128_SHORT_FORMAT "%+.16Le"

inline void printnum(long double rtmp) {
    printf(FLOAT128_FORMAT, rtmp);
    return;
}
inline void printnum(std::complex<long double> ctmp) {
    printf(FLOAT128_FORMAT FLOAT128_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}
inline void sprintnum(char *buf, long double rtmp) {
    snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, rtmp);
    return;
}
inline void sprintnum(char *buf, std::complex<long double> ctmp) {
    snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT FLOAT128_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}
inline void printnum_short(long double rtmp) {
    printf(FLOAT128_SHORT_FORMAT, rtmp);
    return;
}
inline void printnum_short(std::complex<long double> ctmp) {
    printf(FLOAT128_SHORT_FORMAT FLOAT128_SHORT_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}
inline void sprintnum_short(char *buf, long double rtmp) {
    snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, rtmp);
    return;
}
inline void sprintnum_short(char *buf, std::complex<long double> ctmp) {
    snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT FLOAT128_SHORT_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include <math.h>
#include <complex.h>
#include <complex>
#ifdef __cplusplus
}
#endif
typedef long double _Float128;

// implementation of sign transfer function.
inline long double sign(long double a, long double b) {
    long double mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline long double castREAL__Float128(mplapackint n) {
    long double ret = n;
    return ret;
}

inline mplapackint castINTEGER__Float128(long double a) {
    mplapackint i = a;
    return i;
}

inline long nint(long double a) {
    long i;
    long double tmp;
    a = a + 0.5;
    tmp = floorl(a);
    i = (long)tmp;
    return i;
}

inline double cast2double(long double a) { return (double)a; }
inline long double pi(long double dummy) { return M_PIl; }

#endif //___MPLAPACK_LONGDOUBLE_IS_BINARY128___

#if defined ___MPLAPACK_WANT_LIBQUADMATH___

#ifdef __cplusplus
extern "C" {
#endif
#include <quadmath.h>
#ifdef __cplusplus
}
#endif

#if defined(__clang__)
#  if !__has_keyword(_Float128)
     typedef __float128 _Float128;
#  endif
#elif defined(__GNUC__)
/* Skip typedef for GCC 13+ to avoid redeclaration error */
#  if (__GNUC__ < 13)
     typedef __float128 _Float128;
#  endif
#else
   typedef __float128 _Float128;
#endif

#if defined ___MPLAPACK_INTERNAL___
#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif
#include <string.h>

#define FLOAT128_FORMAT "%+-#*.40Qe"
#define FLOAT128_SHORT_FORMAT "%+-#*.16Qe"

inline void printnum(_Float128 rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__];
    int n = quadmath_snprintf(buf, sizeof buf, FLOAT128_FORMAT, width, rtmp);
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    return;
}

inline void printnum(std::complex<_Float128> rtmp) {
    int width = 42, n;
    char buf[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    n = quadmath_snprintf(buf, sizeof buf, FLOAT128_FORMAT, width, rtmp.real());
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    n = quadmath_snprintf(buf2, sizeof buf, FLOAT128_FORMAT, width, rtmp.imag());
    if ((size_t)n < sizeof buf2)
        printf("%s", buf2);
    printf("i");
    return;
}

inline void sprintnum(char *buf, _Float128 rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT, width, rtmp);
    return;
}

inline void sprintnum(char *buf, std::complex<_Float128> rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_FORMAT FLOAT128_FORMAT, width, rtmp.real(), rtmp.imag());
    return;
}

inline void printnum_short(_Float128 rtmp) {
    int width = 42;
    char buf[__MPLAPACK_BUFLEN__];
    int n = quadmath_snprintf(buf, sizeof buf, FLOAT128_SHORT_FORMAT, width, rtmp);
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    return;
}

inline void printnum_short(std::complex<_Float128> rtmp) {
    int width = 42, n;
    char buf[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    n = quadmath_snprintf(buf, sizeof buf, FLOAT128_SHORT_FORMAT, width, rtmp.real());
    if ((size_t)n < sizeof buf)
        printf("%s", buf);
    n = quadmath_snprintf(buf2, sizeof buf, FLOAT128_SHORT_FORMAT, width, rtmp.imag());
    if ((size_t)n < sizeof buf2)
        printf("%s", buf2);
    printf("i");
    return;
}

inline void sprintnum_short(char *buf, _Float128 rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT, width, rtmp);
    return;
}

inline void sprintnum_short(char *buf, std::complex<_Float128> rtmp) {
    int width = 42;
    quadmath_snprintf(buf, __MPLAPACK_BUFLEN__, FLOAT128_SHORT_FORMAT FLOAT128_SHORT_FORMAT, width, rtmp.real(), rtmp.imag());
    return;
}

#endif

inline _Float128 pow(const _Float128 &a, const _Float128 &b) { return powq(a, b); }
inline _Float128 pow(const long &a, const long &b) { return powq((_Float128)a, (_Float128)b); }
inline _Float128 pow(const int &a, const long &b) { return powq((_Float128)a, (_Float128)b); }
inline _Float128 pow(const _Float128 &a, const long &b) { return powq(a, (_Float128)b); }

inline _Float128 sqrt(const _Float128 &a) { return sqrtq(a); }

inline _Float128 sin(_Float128 a) { return sinq(a); }
inline _Float128 sinh(_Float128 a) { return sinhq(a); }
inline _Float128 cos(_Float128 a) { return cosq(a); }
inline _Float128 cosh(_Float128 a) { return coshq(a); }

inline _Float128 atan2(_Float128 a, _Float128 b) { return atan2q(a, b); }

inline _Float128 exp(const _Float128 &a) { return expq(a); }
inline _Float128 log(const _Float128 &a) { return logq(a); }
inline _Float128 log10(const _Float128 &a) { return log10q(a); }
inline _Float128 log2(const _Float128 &a) { return logq(a) / logq(2.0); }
inline _Float128 ldexp(const _Float128 &a, int exp) { return ldexpq(a, exp); }
inline _Float128 nextafter(const _Float128 &a, const _Float128 &b) { return nextafterq(a, b); }

#ifdef __cplusplus
extern "C" {
#endif
#include <complex.h>
#include <complex>
#include <math.h>
#ifdef __cplusplus
}
#endif

inline _Float128 abs(const std::complex<_Float128> &a) {
    _Float128 _Complex b;
    _Float128 c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    c = cabsq(b);
    return c;
}

inline std::complex<_Float128> sqrt(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csqrtq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> sin(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = csinq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> cos(const std::complex<_Float128> a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = ccosq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> exp(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = cexpq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

inline std::complex<_Float128> log(const std::complex<_Float128> &a) {
    _Float128 _Complex b, tmp;
    std::complex<_Float128> c;
    __real__(b) = (a.real());
    __imag__(b) = (a.imag());
    tmp = clogq(b);
    c.real(__real__(tmp));
    c.imag(__imag__(tmp));
    return c;
}

// implementation of sign transfer function.
inline _Float128 sign(_Float128 a, _Float128 b) {
    _Float128 mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline _Float128 castREAL__Float128(mplapackint n) {
    _Float128 ret = n;
    return ret;
}

inline mplapackint castINTEGER__Float128(_Float128 a) {
    mplapackint i = a;
    return i;
}

inline long nint(_Float128 a) {
    long i;
    _Float128 tmp;
    a = a + 0.5;
    tmp = floorq(a);
    i = (long)tmp;
    return i;
}

inline _Float128 ceil(_Float128 a) { return ceilq(a); }
inline double cast2double(_Float128 a) { return (double)a; }

inline _Float128 pi(_Float128 dummy) { return M_PIq; }

#endif //___MPLAPACK_WANT_LIBQUADMATH___

// Following specialization should be done in the libc/compiler side.
#if !defined ___MPLAPACK__FLOAT128_IS_LONGDOUBLE___ && !defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
template <> template <> inline std::complex<_Float128> &std::complex<_Float128>::operator/=(const std::complex<_Float128> &b) {
    _Float128 abr, abi, ratio, den;
    if ((abr = b.real()) < 0.)
        abr = -abr;
    if ((abi = b.imag()) < 0.)
        abi = -abi;
    if (abr <= abi) {
        if (abi == 0) {
            if (this->imag() != 0 || this->real() != 0)
                abi = 1.;
            *this = std::complex<_Float128>(abi / abr, abi / abr);
            return (*this);
        }
        ratio = b.real() / b.imag();
        den = b.imag() * (1.0 + ratio * ratio);
        (*this) = std::complex<_Float128>((this->real() * ratio + this->imag()) / den, (this->imag() * ratio - this->real()) / den);
    } else {
        ratio = b.imag() / b.real();
        den = b.real() * (1.0 + ratio * ratio);
        (*this) = std::complex<_Float128>((this->real() + this->imag() * ratio) / den, (this->imag() - this->real() * ratio) / den);
    }
    return *this;
}
#endif

static inline _Float128 cabs1(const std::complex<_Float128> &z) { return abs(z.real()) + abs(z.imag()); }

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

#ifndef MPLAPACK_MINMAX_FLOAT128_DEFINED
#define MPLAPACK_MINMAX_FLOAT128_DEFINED

inline _Float128 min(_Float128 a, _Float128 b) { return (a > b) ? b : a; }
inline _Float128 max(_Float128 a, _Float128 b) { return (a < b) ? b : a; }

inline _Float128 min(_Float128 a, _Float128 b, _Float128 c) {
    _Float128 r = min(a, b);
    return min(r, c);
}
inline _Float128 max(_Float128 a, _Float128 b, _Float128 c) {
    _Float128 r = max(a, b);
    return max(r, c);
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<_Float128, std::decay_t<Args>> && ...)>> inline _Float128 min(_Float128 a, _Float128 b, _Float128 c, Args... rest) {
    _Float128 r = min(a, b, c);
    ((r = min(r, rest)), ...);
    return r;
}

template <typename... Args, typename = std::enable_if_t<(std::is_same_v<_Float128, std::decay_t<Args>> && ...)>> inline _Float128 max(_Float128 a, _Float128 b, _Float128 c, Args... rest) {
    _Float128 r = max(a, b, c);
    ((r = max(r, rest)), ...);
    return r;
}

#endif // MPLAPACK_MINMAX_FLOAT128_DEFINED

#include <type_traits>

// Square for INTEGER (workspace sizes, indices).
#ifndef MPLAPACK_POW2_MPLAPACKINT_DEFINED
#define MPLAPACK_POW2_MPLAPACKINT_DEFINED
inline mplapackint pow2(mplapackint a) { return a * a; }
#endif // MPLAPACK_POW2_MPLAPACKINT_DEFINED

// Square and quartic for REAL (computational kernels)
#ifndef MPLAPACK_POW2_POW4_FLOAT128_DEFINED
#define MPLAPACK_POW2_POW4_FLOAT128_DEFINED
inline _Float128 pow2(const _Float128 &a) { return a * a; }
inline _Float128 pow4(const _Float128 &a) {
    _Float128 tmp = a * a;
    return tmp * tmp;
}
inline std::complex<_Float128> pow2(const std::complex<_Float128> &a) { return a * a; }
inline std::complex<_Float128> pow4(const std::complex<_Float128> &a) {
    std::complex<_Float128> tmp = a * a;
    return tmp * tmp;
}
#endif // MPLAPACK_POW2_POW4_FLOAT128_DEFINED

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

// Integer ceil for _Float128.
// Returns ceil(x) as mplapackint.
#ifndef MPLAPACK_ICEIL_FLOAT128_DEFINED
#define MPLAPACK_ICEIL_FLOAT128_DEFINED
inline mplapackint iceil(_Float128 x) {
    mplapackint t = static_cast<mplapackint>(x); // trunc toward zero
    if (x > static_cast<_Float128>(t)) {
        ++t;
    }
    return t;
}
#endif // MPLAPACK_ICEIL_FLOAT128_DEFINED

#ifndef MPLAPACK_MOD_UTILS_H
#define MPLAPACK_MOD_UTILS_H
inline mplapackint mod(mplapackint a, mplapackint b) { return a % b; }
#endif // MPLAPACK_MOD_UTILS_H

#endif

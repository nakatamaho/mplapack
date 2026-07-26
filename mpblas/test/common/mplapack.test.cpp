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

#define __MPLAPACK_BUILD_DEBUG_CPP__

#include <mpblas.h>
#include <mplapack_compare_debug.h>

void __attribute__((constructor)) mplapack_debug_initialize(void);
void mplapack_debug_initialize(void) { gmp_randinit_default(uniformrandomstate_mpfr); }

void __attribute__((destructor)) mplapack_debug_finalize(void);
void mplapack_debug_finalize(void) {
    gmp_randclear(uniformrandomstate_mpfr);
    mpfr_free_cache(); // we always use MPFR when debugging
}

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void __attribute__((constructor)) mplapack_debug_initialize_gmp(void);
void mplapack_debug_initialize_gmp(void) { uniformrandomstate_gmp = new gmp_randclass(gmp_randinit_default); }
#endif

mpfr_class mpf_randomnumber(mpfr_class dummy) {
    mpfr_class mtmp;

    mtmp = urandomb(uniformrandomstate_mpfr);
    mtmp = 2.0 * mtmp - 1.0;

    return mtmp;
}

mpc_class mpc_randomnumber(mpc_class dummy) {
    mpc_class ctmp = urandom_c(uniformrandomstate_mpfr);
    return ctmp;
}

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
double mpf_randomnumber(double dummy) {
    double mtmp = drand48();
    return mtmp;
}

complex<double> mpc_randomnumber(complex<double>) {
    std::complex<double> ctmp;
    ctmp.real(drand48());
    ctmp.imag(drand48());
    return ctmp;
}

void set_random_number(double &a, mpfr_class &b) {
    double dummy = 0.0;
    a = mpf_randomnumber(dummy);
    b = a;
}

void set_random_number(complex<double> &a, mpc_class &b) {
    complex<double> dummy, p;
    a = mpc_randomnumber(dummy);
    b = a;
}

void set_random_number1to2(double &a, mpfr_class &b) {
    double dummy = 0.0;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = a;
}

void set_random_number1to2(complex<double> &a, mpc_class &b) {
    complex<double> dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0)
        p = 1.0;
    else
        p = -1.0;
    if (a.imag() > 0.0)
        q = 1.0;
    else
        q = -1.0;
    a = a + complex<double>(p, q);
    b = a;
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
mpf_class mpf_randomnumber(mpf_class dummy) {
    mpf_class mtmp;

    mtmp = uniformrandomstate_gmp->get_f();
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}

mpfc_class mpc_randomnumber(mpfc_class dummy) {
    mpfc_class ctmp;
    mpf_class mtmp1;
    mpf_class mtmp2;

    mtmp1 = uniformrandomstate_gmp->get_f();
    mtmp1 = 2.0 * mtmp1 - 1.0;

    mtmp2 = uniformrandomstate_gmp->get_f();
    mtmp2 = 2.0 * mtmp2 - 1.0;

    ctmp.real(mtmp1);
    ctmp.imag(mtmp2);
    return ctmp;
}

void set_random_number(mpfr_class &a, mpf_class &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    b = cast2mpf_class(a);
}

void set_random_number(mpc_class &a, mpfc_class &b) {
    mpc_class dummy;
    a = mpc_randomnumber(dummy);
    b = cast2mpc_class(a);
}

void set_random_number1to2(mpfr_class &a, mpf_class &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = cast2mpf_class(a);
}

void set_random_number1to2(mpc_class &a, mpfc_class &b) {
    mpc_class dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0)
        p = 1.0;
    else
        p = -1.0;
    if (a.imag() > 0.0)
        q = 1.0;
    else
        q = -1.0;
    a = a + complex<double>(p, q);
    b = cast2mpc_class(a);
}

#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
qd_real mpf_randomnumber(qd_real dummy) {
    qd_real mtmp;

    mtmp = qdrand(); // uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}

qd_complex mpc_randomnumber(qd_complex dummy) {
    qd_complex ctmp;
    qd_real mtmp1;
    qd_real mtmp2;

    mtmp1 = qdrand();
    mtmp1 = 2.0 * mtmp1 - 1.0;

    mtmp2 = qdrand();
    mtmp2 = 2.0 * mtmp2 - 1.0;

    ctmp.real() = mtmp1;
    ctmp.imag() = mtmp2;
    return ctmp;
}

void set_random_number(mpfr_class &a, qd_real &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    b = cast2qd_real(a);
}

void set_random_number(mpc_class &a, qd_complex &b) {
    mpc_class dummy;
    a = mpc_randomnumber(dummy);
    b = cast2qd_complex(a);
}

void set_random_number1to2(mpfr_class &a, qd_real &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = cast2qd_real(a);
}

void set_random_number1to2(mpc_class &a, qd_complex &b) {
    mpc_class dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0)
        p = 1.0;
    else
        p = -1.0;
    if (a.imag() > 0.0)
        q = 1.0;
    else
        q = -1.0;
    a = a + complex<double>(p, q);
    b = cast2qd_complex(a);
}

#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
dd_real mpf_randomnumber(dd_real dummy) {
    dd_real mtmp;

    mtmp = ddrand(); // uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}

dd_complex mpc_randomnumber(dd_complex dummy) {
    dd_complex ctmp;
    dd_real mtmp1;
    dd_real mtmp2;

    mtmp1 = ddrand();
    mtmp1 = 2.0 * mtmp1 - 1.0;

    mtmp2 = ddrand();
    mtmp2 = 2.0 * mtmp2 - 1.0;

    ctmp.real() = mtmp1;
    ctmp.imag() = mtmp2;
    return ctmp;
}

void set_random_number(mpfr_class &a, dd_real &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    b = cast2dd_real(a);
}

void set_random_number(mpc_class &a, dd_complex &b) {
    mpc_class dummy;
    a = mpc_randomnumber(dummy);
    b = cast2dd_complex(a);
}

void set_random_number1to2(mpfr_class &a, dd_real &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = cast2dd_real(a);
}

void set_random_number1to2(mpc_class &a, dd_complex &b) {
    mpc_class dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0)
        p = 1.0;
    else
        p = -1.0;
    if (a.imag() > 0.0)
        q = 1.0;
    else
        q = -1.0;
    a = a + complex<double>(p, q);
    b = cast2dd_complex(a);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
double mpf_randomnumber(double dummy) {
    double mtmp = drand48();
    return mtmp;
}

complex<double> mpc_randomnumber(complex<double>) {
    std::complex<double> ctmp;
    ctmp.real(drand48());
    ctmp.imag(drand48());
    return ctmp;
}

void set_random_number(mpfr_class &a, double &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    b = a;
}

void set_random_number(mpc_class &a, complex<double> &b) {
    mpc_class dummy;
    a = mpc_randomnumber(dummy);
    b = a;
}

void set_random_number1to2(mpfr_class &a, double &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = a;
}

void set_random_number1to2(mpc_class &a, complex<double> &b) {
    mpc_class dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0)
        p = 1.0;
    else
        p = -1.0;
    if (a.imag() > 0.0)
        q = 1.0;
    else
        q = -1.0;
    a = a + complex<double>(p, q);
    b = a;
}

#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
#include <random>
#include <complex>
#include <cstdint>
#include <cmath>
#include <cfloat>

// Optional sanity checks (fail fast on misconfig).
#if defined(MPLAPACK_BINARY80_MODE) && defined(MPLAPACK_BINARY80_MODE_LDBL80)
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
static_assert(LDBL_MANT_DIG == 64, "long double is not IEEE binary80 (mantissa must be 64 bits).");
#endif
#endif

#if defined(MPLAPACK_BINARY80_MODE) && defined(MPLAPACK_BINARY80_MODE_FLOAT64X)
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
static_assert(FLT64X_MANT_DIG == 64, "_Float64x is not IEEE binary80 (mantissa must be 64 bits).");
#endif
#endif

static inline std::mt19937_64 &mplapack_rng64_binary80() {
    // Deterministic by default; change seed policy if needed.
    static thread_local std::mt19937_64 rng(0x9e3779b97f4a7c15ULL);
    return rng;
}

// Uniform in [0, 1) with 64-bit resolution (binary80 significand bits).
static inline mplapack_binary80_t uniform01_binary80() {
    constexpr int kFracBits = 64;

    uint64_t r = mplapack_rng64_binary80()(); // 64 random bits
    // Scale by 2^-64 using a power-of-two factor (exact).
    const long double scale_ld = std::ldexp(1.0L, -kFracBits);
    const mplapack_binary80_t scale = static_cast<mplapack_binary80_t>(scale_ld);

    // r is in [0, 2^64-1], representable exactly in binary80.
    return static_cast<mplapack_binary80_t>(r) * scale; // [0,1)
}

mplapack_binary80_t mpf_randomnumber(mplapack_binary80_t /*dummy*/) {
    const mplapack_binary80_t u = uniform01_binary80();                                   // [0,1)
    return static_cast<mplapack_binary80_t>(2) * u - static_cast<mplapack_binary80_t>(1); // [-1,1)
}

std::complex<mplapack_binary80_t> mpc_randomnumber(std::complex<mplapack_binary80_t> /*dummy*/) {
    const mplapack_binary80_t re = mpf_randomnumber(static_cast<mplapack_binary80_t>(0));
    const mplapack_binary80_t im = mpf_randomnumber(static_cast<mplapack_binary80_t>(0));
    return std::complex<mplapack_binary80_t>(re, im);
}

void set_random_number(mpfr_class &a, mplapack_binary80_t &b) {
    a = mpf_randomnumber(static_cast<mplapack_binary80_t>(0));
    b = cast2binary80_t(a);
}

void set_random_number(mpc_class &a, std::complex<mplapack_binary80_t> &b) {
    a = mpc_randomnumber(std::complex<mplapack_binary80_t>(0, 0));
    b.real(cast2binary80_t(a.real()));
    b.imag(cast2binary80_t(a.imag()));
}

void set_random_number1to2(mpfr_class &a, mplapack_binary80_t &b) {
    a = mpf_randomnumber(static_cast<mplapack_binary80_t>(0));
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = cast2binary80_t(a);
}

void set_random_number1to2(mpc_class &a, std::complex<mplapack_binary80_t> &b) {
    a = mpc_randomnumber(std::complex<mplapack_binary80_t>(0, 0));

    const mpfr_class p = (a.real() > 0.0) ? mpfr_class(1.0) : mpfr_class(-1.0);
    const mpfr_class q = (a.imag() > 0.0) ? mpfr_class(1.0) : mpfr_class(-1.0);

    a = a + mpc_class(p, q);

    b.real(cast2binary80_t(a.real()));
    b.imag(cast2binary80_t(a.imag()));
}
#endif //___MPLAPACK_BUILD_WITH_BINARY80___

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
#include <random>
static inline std::mt19937_64 &mplapack_rng64() {
    // Deterministic seed; change if you want different streams.
    static thread_local std::mt19937_64 rng(0x243f6a8885a308d3ULL);
    return rng;
}
// Uniform in [0, 1) with 113-bit resolution (binary128 significand bits).
static inline __float128 uniform01_binary128() {
    // Generates a uniform random value in [0, 1) as an integer multiple of 2^-113,
    // using the full 113-bit significand of __float128.
    // Distribution: integer-grid uniform (not float-space uniform).
    constexpr int kFracBits = 113;
    constexpr int lo_bits = kFracBits - 64; // 49

    auto &rng = mplapack_rng64(); // ensure single RNG instance
    uint64_t hi = rng();
    uint64_t lo_top = rng() >> (64 - lo_bits); // top 49 bits

    __float128 r = static_cast<__float128>(hi);
    r = ldexpq(r, lo_bits); // r = hi * 2^49
    r += static_cast<__float128>(lo_top);
    r = ldexpq(r, -kFracBits); // r /= 2^113
    return r;
}
mplapack_binary128_t mpf_randomnumber(mplapack_binary128_t /*dummy*/) {
    __float128 u = uniform01_binary128();                      // [0,1)
    return static_cast<mplapack_binary128_t>(2.0Q * u - 1.0Q); // [-1,1)
}
std::complex<mplapack_binary128_t> mpc_randomnumber(std::complex<mplapack_binary128_t> /*dummy*/) {
    mplapack_binary128_t re = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    mplapack_binary128_t im = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    return std::complex<mplapack_binary128_t>(re, im);
}
void set_random_number(mpfr_class &a, mplapack_binary128_t &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    b = cast2binary128_t(a);
}
void set_random_number(mpc_class &a, std::complex<mplapack_binary128_t> &b) {
    mpc_class dummy;
    a = mpc_randomnumber(dummy);
    b.real(cast2binary128_t(a.real()));
    b.imag(cast2binary128_t(a.imag()));
}
void set_random_number1to2(mpfr_class &a, mplapack_binary128_t &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = cast2binary128_t(a);
}
void set_random_number1to2(mpc_class &a, std::complex<mplapack_binary128_t> &b) {
    mpc_class dummy;
    mplapack_binary128_t p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0)
        p = 1.0;
    else
        p = -1.0;
    if (a.imag() > 0.0)
        q = 1.0;
    else
        q = -1.0;
    a = a + std::complex<mplapack_binary128_t>(p, q);
    b.real(cast2binary128_t(a.real()));
    b.imag(cast2binary128_t(a.imag()));
}
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
#include <random>
#include <complex>
#include <cstdint>
#include <math.h> // ldexpf128 / scalbnf128 (glibc)

static inline std::mt19937_64 &mplapack_rng64() {
    // Deterministic seed; change if you want different streams.
    static thread_local std::mt19937_64 rng(0x243f6a8885a308d3ULL);
    return rng;
}

// Prefer ldexpf128; fall back to scalbnf128 if needed.
static inline _Float128 mplapack_ldexp_f128(_Float128 x, int e) {
#if defined(ldexpf128) || (defined(__GLIBC__) && !defined(__cplusplus))
    // Some libcs expose it as a macro; keep it simple.
    return ::ldexpf128(x, e);
#else
    // Many libm provide scalbnf128 even if ldexpf128 is missing.
    return ::scalbnf128(x, e);
#endif
}

// Uniform in [0, 1) with 113-bit resolution (binary128 significand bits).
static inline _Float128 uniform01_binary128_f128() {
    constexpr int kFracBits = 113;
    constexpr int kShift = 128 - kFracBits; // 15

    uint64_t hi = mplapack_rng64()();
    uint64_t lo = mplapack_rng64()();

    // Extract top 113 bits of (hi:lo) by right-shifting the 128-bit value by 15.
    uint64_t a = hi >> kShift;                           // top 49 bits
    uint64_t b = (hi << (64 - kShift)) | (lo >> kShift); // next 64 bits

    // Build exact integer: x = a * 2^64 + b  (0 <= x < 2^113)
    _Float128 r = mplapack_ldexp_f128(static_cast<_Float128>(a), 64); // a * 2^64
    r += static_cast<_Float128>(b);                                   // + b

    // Scale to [0, 1)
    r = mplapack_ldexp_f128(r, -kFracBits); // / 2^113
    return r;
}

mplapack_binary128_t mpf_randomnumber(mplapack_binary128_t /*dummy*/) {
    _Float128 u = uniform01_binary128_f128();            // [0,1)
    _Float128 v = (_Float128)2.0L * u - (_Float128)1.0L; // [-1,1)
    return static_cast<mplapack_binary128_t>(v);
}

std::complex<mplapack_binary128_t> mpc_randomnumber(std::complex<mplapack_binary128_t> /*dummy*/) {
    mplapack_binary128_t re = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    mplapack_binary128_t im = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    return std::complex<mplapack_binary128_t>(re, im);
}

void set_random_number(mpfr_class &a, mplapack_binary128_t &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    b = cast2binary128_t(a);
}

void set_random_number(mpc_class &a, std::complex<mplapack_binary128_t> &b) {
    mpc_class dummy;
    a = mpc_randomnumber(dummy);
    b.real(cast2binary128_t(a.real()));
    b.imag(cast2binary128_t(a.imag()));
}

void set_random_number1to2(mpfr_class &a, mplapack_binary128_t &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = cast2binary128_t(a);
}

void set_random_number1to2(mpc_class &a, std::complex<mplapack_binary128_t> &b) {
    mpc_class dummy;
    mplapack_binary128_t p, q;
    a = mpc_randomnumber(dummy);

    p = (a.real() > 0.0) ? (mplapack_binary128_t)1.0L : (mplapack_binary128_t)(-1.0L);
    q = (a.imag() > 0.0) ? (mplapack_binary128_t)1.0L : (mplapack_binary128_t)(-1.0L);

    a = a + std::complex<mplapack_binary128_t>(p, q);
    b.real(cast2binary128_t(a.real()));
    b.imag(cast2binary128_t(a.imag()));
}
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
#include <random>
#include <complex>
#include <cstdint>
#include <cmath>
#include <cfloat>

static inline std::mt19937_64 &mplapack_rng64_binary128_ldbl() {
    static thread_local std::mt19937_64 rng(0x243f6a8885a308d3ULL);
    return rng;
}

// Uniform in [0, 1) with LDBL_MANT_DIG-bit resolution.
static inline long double uniform01_binary128_ldbl() {
    constexpr int kFracBits = LDBL_MANT_DIG; // 64 (binary80), 113 (binary128 long double)
    static_assert(kFracBits >= 2 && kFracBits <= 113, "Unexpected LDBL_MANT_DIG");

    uint64_t hi = mplapack_rng64_binary128_ldbl()();
    uint64_t lo = mplapack_rng64_binary128_ldbl()();

    if constexpr (kFracBits <= 64) {
        // Top k bits of the 128-bit stream live entirely in hi.
        // Equivalent to: ((hi<<64 | lo) >> (128-k)) but without 128-bit ints.
        constexpr int sh = 64 - kFracBits;        // 0..63
        uint64_t x = (sh == 0) ? hi : (hi >> sh); // k-bit integer (< 2^k)
        long double r = static_cast<long double>(x);
        return std::ldexp(r, -kFracBits); // r / 2^k
    } else {
        // k in (64,113], so shift s = 128-k is in [15,63].
        constexpr int s = 128 - kFracBits; // 1..63 in practice here

        // 128-bit right shift by s using two 64-bit limbs (hi:lo).
        uint64_t y_hi = hi >> s;
        uint64_t y_lo = (hi << (64 - s)) | (lo >> s);

        // Build exact integer y = y_hi*2^64 + y_lo (needs >=113-bit long double to be exact).
        long double r = std::ldexp(static_cast<long double>(y_hi), 64);
        r += static_cast<long double>(y_lo);
        return std::ldexp(r, -kFracBits); // / 2^k
    }
}

mplapack_binary128_t mpf_randomnumber(mplapack_binary128_t /*dummy*/) {
    long double u = uniform01_binary128_ldbl();                // [0,1)
    return static_cast<mplapack_binary128_t>(2.0L * u - 1.0L); // [-1,1)
}

std::complex<mplapack_binary128_t> mpc_randomnumber(std::complex<mplapack_binary128_t> /*dummy*/) {
    mplapack_binary128_t re = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    mplapack_binary128_t im = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    return std::complex<mplapack_binary128_t>(re, im);
}

void set_random_number(mpfr_class &a, mplapack_binary128_t &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    b = cast2binary128_t(a);
}

void set_random_number(mpc_class &a, std::complex<mplapack_binary128_t> &b) {
    mpc_class dummy;
    a = mpc_randomnumber(dummy);
    b.real(cast2binary128_t(a.real()));
    b.imag(cast2binary128_t(a.imag()));
}

void set_random_number1to2(mpfr_class &a, mplapack_binary128_t &b) {
    mpfr_class dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0)
        a = a + 1.0;
    else
        a = a - 1.0;
    b = cast2binary128_t(a);
}

void set_random_number1to2(mpc_class &a, std::complex<mplapack_binary128_t> &b) {
    mpc_class dummy;
    mplapack_binary128_t p, q;
    a = mpc_randomnumber(dummy);
    p = (a.real() > 0.0) ? (mplapack_binary128_t)1.0L : (mplapack_binary128_t)(-1.0L);
    q = (a.imag() > 0.0) ? (mplapack_binary128_t)1.0L : (mplapack_binary128_t)(-1.0L);
    a = a + std::complex<mplapack_binary128_t>(p, q);
    b.real(cast2binary128_t(a.real()));
    b.imag(cast2binary128_t(a.imag()));
}
#endif // MPLAPACK_BINARY128_MODE
#endif // ___MPLAPACK_BUILD_WITH_BINARY128___

REAL_REF infnorm(COMPLEX_REF *vec_ref, COMPLEX *vec, int len, int inc) {
    COMPLEX_REF ctemp;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
        ctemp = vec[i] - vec_ref[i];
        if (inorm < abs(ctemp)) {
            inorm = abs(ctemp);
        }
    }
    return inorm;
}

REAL_REF infnorm(REAL_REF *vec_ref, REAL *vec, int len, int inc) {
    REAL_REF ctemp;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
        ctemp = vec[i] - vec_ref[i];
        if (inorm < abs(ctemp)) {
            inorm = abs(ctemp);
        }
    }
    return inorm;
}

REAL_REF infnorm_mat(REAL_REF *vec_ref, REAL *vec, int n, int m, int lda) {
    REAL_REF ctemp;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            ctemp = vec[i + j * lda] - vec_ref[i + j * lda];
            if (inorm < abs(ctemp)) {
                inorm = abs(ctemp);
            }
        }
    }
    return inorm;
}

REAL_REF reldiff(REAL_REF *vec_ref, REAL *vec, int len, int inc) {
    REAL_REF ctemp1, ctemp2, ctemp3;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
        ctemp1 = (REAL_REF)abs(vec[i]);
        ctemp2 = abs(vec_ref[i]);
        ctemp3 = max(ctemp1, ctemp2);
        ctemp1 = vec[i] - vec_ref[i];
        if (inorm < abs(ctemp1) / ctemp3) {
            inorm = abs(ctemp1) / ctemp3;
        }
    }
    return inorm;
}

INTEGER_REF infnorm(INTEGER_REF *vec_ref, INTEGER *vec, int len, int inc) {
    INTEGER_REF inorm = 0;
    INTEGER_REF temp = 0;
    for (int i = 0; i < len * inc; i = i + inc) {
        temp = vec[i] - vec_ref[i];
        if (inorm < abs(temp)) {
            inorm = abs(temp);
        }
    }
    return inorm;
}

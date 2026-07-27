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

void __attribute__((destructor)) mplapack_debug_finalize(void);
void mplapack_debug_finalize(void) {
    mpfr_free_cache(); // we always use MPFR when debugging
}

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
mpfr_class mpf_randomnumber(mpfr_class /*dummy*/) {
    mpfr_class mtmp;

    mtmp = uniformrandomstate_mpfr.get_fr();
    mtmp = 2.0 * mtmp - 1.0;

    return mtmp;
}

mpc_class mpc_randomnumber(mpc_class /*dummy*/) {
    mpfr_class real_part;
    mpfr_class imag_part;
    real_part = uniformrandomstate_mpfr.get_fr();
    imag_part = uniformrandomstate_mpfr.get_fr();
    real_part = 2.0 * real_part - 1.0;
    imag_part = 2.0 * imag_part - 1.0;
    return mpc_class(real_part, imag_part);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
mpf_class mpf_randomnumber(mpf_class /*dummy*/) {
    mpf_class value;
    value = uniformrandomstate_gmp.get_f();
    return 2.0 * value - 1.0;
}

mpfc_class mpc_randomnumber(mpfc_class /*dummy*/) {
    mpf_class real_part;
    mpf_class imag_part;
    real_part = uniformrandomstate_gmp.get_f();
    imag_part = uniformrandomstate_gmp.get_f();
    return mpfc_class(2.0 * real_part - 1.0, 2.0 * imag_part - 1.0);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
qd_real mpf_randomnumber(qd_real /*dummy*/) {
    return 2.0 * qdrand() - 1.0;
}

qd_complex mpc_randomnumber(qd_complex /*dummy*/) {
    return qd_complex(2.0 * qdrand() - 1.0, 2.0 * qdrand() - 1.0);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
dd_real mpf_randomnumber(dd_real /*dummy*/) {
    return 2.0 * ddrand() - 1.0;
}

dd_complex mpc_randomnumber(dd_complex /*dummy*/) {
    return dd_complex(2.0 * ddrand() - 1.0, 2.0 * ddrand() - 1.0);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
double mpf_randomnumber(double /*dummy*/) {
    return 2.0 * drand48() - 1.0;
}

complex<double> mpc_randomnumber(complex<double> /*dummy*/) {
    return complex<double>(mpf_randomnumber(0.0), mpf_randomnumber(0.0));
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdint>
#include <random>

#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
static_assert(LDBL_MANT_DIG == 64, "long double is not IEEE binary80");
#elif MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
static_assert(FLT64X_MANT_DIG == 64, "_Float64x is not IEEE binary80");
#endif

static std::mt19937_64 &mplapack_rng64_binary80() {
    static thread_local std::mt19937_64 rng(0x9e3779b97f4a7c15ULL);
    return rng;
}

static mplapack_binary80_t uniform01_binary80() {
    constexpr int fraction_bits = 64;
    const uint64_t bits = mplapack_rng64_binary80()();
    const long double scale = std::ldexp(1.0L, -fraction_bits);
    return static_cast<mplapack_binary80_t>(bits) * static_cast<mplapack_binary80_t>(scale);
}

mplapack_binary80_t mpf_randomnumber(mplapack_binary80_t /*dummy*/) {
    const mplapack_binary80_t value = uniform01_binary80();
    return static_cast<mplapack_binary80_t>(2) * value - static_cast<mplapack_binary80_t>(1);
}

complex<mplapack_binary80_t> mpc_randomnumber(complex<mplapack_binary80_t> /*dummy*/) {
    const mplapack_binary80_t real_part = mpf_randomnumber(static_cast<mplapack_binary80_t>(0));
    const mplapack_binary80_t imag_part = mpf_randomnumber(static_cast<mplapack_binary80_t>(0));
    return complex<mplapack_binary80_t>(real_part, imag_part);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdint>
#include <random>

static std::mt19937_64 &mplapack_rng64_binary128() {
    static thread_local std::mt19937_64 rng(0x243f6a8885a308d3ULL);
    return rng;
}

#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
static mplapack_binary128_t uniform01_binary128() {
    constexpr int fraction_bits = 113;
    constexpr int low_bits = fraction_bits - 64;
    const uint64_t high = mplapack_rng64_binary128()();
    const uint64_t low = mplapack_rng64_binary128()() >> (64 - low_bits);
    __float128 value = ldexpq(static_cast<__float128>(high), low_bits);
    value += static_cast<__float128>(low);
    return static_cast<mplapack_binary128_t>(ldexpq(value, -fraction_bits));
}
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
static _Float128 mplapack_ldexp_float128(_Float128 value, int exponent) {
#if defined(ldexpf128) || (defined(__GLIBC__) && !defined(__cplusplus))
    return ::ldexpf128(value, exponent);
#else
    return ::scalbnf128(value, exponent);
#endif
}

static mplapack_binary128_t uniform01_binary128() {
    constexpr int fraction_bits = 113;
    constexpr int shift = 128 - fraction_bits;
    const uint64_t high = mplapack_rng64_binary128()();
    const uint64_t low = mplapack_rng64_binary128()();
    const uint64_t upper = high >> shift;
    const uint64_t lower = (high << (64 - shift)) | (low >> shift);
    _Float128 value = mplapack_ldexp_float128(static_cast<_Float128>(upper), 64);
    value += static_cast<_Float128>(lower);
    return static_cast<mplapack_binary128_t>(mplapack_ldexp_float128(value, -fraction_bits));
}
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
static mplapack_binary128_t uniform01_binary128() {
    constexpr int fraction_bits = LDBL_MANT_DIG;
    static_assert(fraction_bits >= 2 && fraction_bits <= 113, "unexpected LDBL_MANT_DIG");
    const uint64_t high = mplapack_rng64_binary128()();
    const uint64_t low = mplapack_rng64_binary128()();
    if constexpr (fraction_bits <= 64) {
        constexpr int shift = 64 - fraction_bits;
        const uint64_t bits = (shift == 0) ? high : (high >> shift);
        return static_cast<mplapack_binary128_t>(std::ldexp(static_cast<long double>(bits), -fraction_bits));
    } else {
        constexpr int shift = 128 - fraction_bits;
        const uint64_t upper = high >> shift;
        const uint64_t lower = (high << (64 - shift)) | (low >> shift);
        long double value = std::ldexp(static_cast<long double>(upper), 64);
        value += static_cast<long double>(lower);
        return static_cast<mplapack_binary128_t>(std::ldexp(value, -fraction_bits));
    }
}
#endif

mplapack_binary128_t mpf_randomnumber(mplapack_binary128_t /*dummy*/) {
    const mplapack_binary128_t value = uniform01_binary128();
    return static_cast<mplapack_binary128_t>(2) * value - static_cast<mplapack_binary128_t>(1);
}

complex<mplapack_binary128_t> mpc_randomnumber(complex<mplapack_binary128_t> /*dummy*/) {
    const mplapack_binary128_t real_part = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    const mplapack_binary128_t imag_part = mpf_randomnumber(static_cast<mplapack_binary128_t>(0));
    return complex<mplapack_binary128_t>(real_part, imag_part);
}
#endif

void set_random_number(REAL_REF &reference_value, REAL &backend_value) {
    backend_value = mpf_randomnumber(backend_value);
    reference_value = cast2ref(backend_value);
}

void set_random_number(COMPLEX_REF &reference_value, COMPLEX &backend_value) {
    backend_value = mpc_randomnumber(backend_value);
    reference_value = cast2ref(backend_value);
}

void set_random_number1to2(REAL_REF &reference_value, REAL &backend_value) {
    backend_value = mpf_randomnumber(backend_value);
    backend_value = backend_value + ((backend_value > REAL(0.0)) ? REAL(1.0) : REAL(-1.0));
    reference_value = cast2ref(backend_value);
}

void set_random_number1to2(COMPLEX_REF &reference_value, COMPLEX &backend_value) {
    backend_value = mpc_randomnumber(backend_value);
    const REAL real_shift = (backend_value.real() > REAL(0.0)) ? REAL(1.0) : REAL(-1.0);
    const REAL imag_shift = (backend_value.imag() > REAL(0.0)) ? REAL(1.0) : REAL(-1.0);
    backend_value = backend_value + COMPLEX(real_shift, imag_shift);
    reference_value = cast2ref(backend_value);
}

REAL_REF infnorm(COMPLEX_REF *vec_ref, COMPLEX *vec, int len, int inc) {
    COMPLEX_REF difference;
    REAL_REF norm = 0.0;
    for (int i = 0; i < len * inc; i += inc) {
        difference = cast2ref(vec[i]) - vec_ref[i];
        if (norm < abs(difference)) {
            norm = abs(difference);
        }
    }
    return norm;
}

REAL_REF infnorm(REAL_REF *vec_ref, REAL *vec, int len, int inc) {
    REAL_REF ctemp;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
        ctemp = cast2ref(vec[i]) - vec_ref[i];
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
            ctemp = cast2ref(vec[i + j * lda]) - vec_ref[i + j * lda];
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
        ctemp1 = abs(cast2ref(vec[i]));
        ctemp2 = abs(vec_ref[i]);
        ctemp3 = max(ctemp1, ctemp2);
        ctemp1 = cast2ref(vec[i]) - vec_ref[i];
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

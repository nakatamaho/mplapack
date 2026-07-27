/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_debug.h,v 1.55 2010/08/18 10:28:48 nakatamaho Exp $
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

#ifndef _MPLAPACK_COMPARE_DEBUG_H_
#define _MPLAPACK_COMPARE_DEBUG_H_

#ifndef ___MPLAPACK_DEBUG_BUILD___
#define ___MPLAPACK_DEBUG_BUILD___
#endif

extern int blas_errno; // for BLAS xerbla dispatch

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using std::complex;
using std::cout;
using std::endl;
using std::max;
using std::min;

#include <mplapack_gmpfrxx_mkII_config.h>
#include <mpfrxx_mkII.h>
#include <mpcxx_mkII.h>
using namespace mpfrxx;

#include <mplapack_print_double.h>

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___
#include <mplapack_gmpfrxx_binary_adapters.h>
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
#define EPSILON 1e-31
#define EPSILON2 1e-27
#define EPSILON3 1e-27
#define EPSILON4 1e-27
#define EPSILON6 1e-24
#define EPSILON7 1e-24
#define EPSILON8 1e-23
#define EPSILON10 1e-22
#define EPSILON11 1e-21
#define EPSILON12 1e-20
#elif defined ___MPLAPACK_BUILD_WITH_MPFR___
#define EPSILON 1e-12
#define EPSILON2 1e-10
#define EPSILON3 1e-9
#define EPSILON4 1e-8
#define EPSILON6 1e-8
#define EPSILON7 1e-7
#define EPSILON8 1e-7
#define EPSILON10 1e-7
#define EPSILON11 1e-7
#define EPSILON12 1e-5
#elif defined ___MPLAPACK_BUILD_WITH_GMP___
#define EPSILON 1e-150
#define EPSILON2 1e-148
#define EPSILON3 1e-148
#define EPSILON4 1e-148
#define EPSILON6 1e-145
#define EPSILON7 1e-145
#define EPSILON8 1e-145
#define EPSILON10 1e-140
#define EPSILON11 1e-140
#define EPSILON12 1e-140
#define EPSILON100 1e-13
#elif defined ___MPLAPACK_BUILD_WITH_QD___
#define EPSILON 1e-61
#define EPSILON2 1e-55
#define EPSILON3 1e-55
#define EPSILON4 1e-55
#define EPSILON6 1e-54
#define EPSILON7 1e-54
#define EPSILON8 1e-54
#define EPSILON10 1e-54
#define EPSILON11 1e-54
#define EPSILON12 1e-54
#elif defined ___MPLAPACK_BUILD_WITH_DD___
#define EPSILON 1e-29
#define EPSILON2 1e-26
#define EPSILON3 1e-26
#define EPSILON4 1e-25
#define EPSILON6 1e-24
#define EPSILON7 1e-23
#define EPSILON8 1e-23
#define EPSILON10 1e-22
#define EPSILON11 1e-21
#define EPSILON12 1e-20
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
#define EPSILON 1e-12
#define EPSILON2 1e-10
#define EPSILON3 1e-9
#define EPSILON4 1e-8
#define EPSILON6 1e-7
#define EPSILON7 1e-7
#define EPSILON8 1e-7
#define EPSILON10 1e-7
#define EPSILON11 1e-7
#define EPSILON12 1e-6
#elif defined ___MPLAPACK_BUILD_WITH_BINARY80___
#define EPSILON 1e-13
#define EPSILON2 1e-11
#define EPSILON3 1e-10
#define EPSILON4 1e-9
#define EPSILON6 1e-8
#define EPSILON7 1e-8
#define EPSILON8 1e-8
#define EPSILON10 1e-8
#define EPSILON11 1e-8
#define EPSILON12 1e-7
#endif

#define TRUE 1
#define FALSE 0

#if defined _WIN32
// no drand48 on mingw
inline double drand48() {
    static const double m_const = 4.6566128730773926e-10; /* = 2^{-31} */
    double d;
    d = std::rand() * m_const;
    return d;
}
#endif

extern int mplapack_errno; // Mxerbla.override.cpp

#if defined __MPLAPACK_BUILD_DEBUG_CPP__
#define _MPLAPACK_DEBUG_EXTERN_
#else
#define _MPLAPACK_DEBUG_EXTERN_ extern
#endif

#if defined(___MPLAPACK_BUILD_WITH_GMP___) || defined(___MPLAPACK_BUILD_WITH_QD___) || defined(___MPLAPACK_BUILD_WITH_DD___) || defined(___MPLAPACK_BUILD_WITH_DOUBLE___) || defined(___MPLAPACK_BUILD_WITH_BINARY80___) || defined(___MPLAPACK_BUILD_WITH_BINARY128___)
#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>
#endif

_MPLAPACK_DEBUG_EXTERN_ mpfrxx::mpfr_randclass uniformrandomstate_mpfr;

#if defined ___MPLAPACK_BUILD_WITH_GMP___
_MPLAPACK_DEBUG_EXTERN_ gmpxx::gmp_randclass uniformrandomstate_gmp;
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
typedef int INTEGER_REF;
typedef double REAL_REF;
typedef complex<double> COMPLEX_REF;
#else
typedef mplapackint INTEGER_REF;
typedef mpfr_class REAL_REF;
typedef mpc_class COMPLEX_REF;
#endif

template <typename T> inline const T &cast2ref(const T &value) { return value; }

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
inline REAL_REF cast2ref(const mpfrxx::mpfr_class &value) { return value.get_d(); }
inline COMPLEX_REF cast2ref(const mpfrxx::mpc_class &value) {
    return COMPLEX_REF(value.real_get_d(), value.imag_get_d());
}
#elif defined ___MPLAPACK_BUILD_WITH_GMP___
inline REAL_REF cast2ref(const gmpxx::mpf_class &value) {
    REAL_REF result = REAL_REF::with_precision(static_cast<mpfr_prec_t>(value.precision()));
    mpfr_set_f(result.mpfr_data(), value.mpf_data(), REAL_REF::default_rounding());
    return result;
}
inline COMPLEX_REF cast2ref(const gmpxx::mpfc_class &value) {
    return COMPLEX_REF(cast2ref(value.real()), cast2ref(value.imag()));
}
#elif defined ___MPLAPACK_BUILD_WITH_QD___
inline REAL_REF cast2ref(const qd_real &value) {
    REAL_REF result = REAL_REF::with_precision(REAL_REF::default_precision());
    mpfr_set_d(result.mpfr_data(), value.x[0], REAL_REF::default_rounding());
    mpfr_add_d(result.mpfr_data(), result.mpfr_data(), value.x[1], REAL_REF::default_rounding());
    mpfr_add_d(result.mpfr_data(), result.mpfr_data(), value.x[2], REAL_REF::default_rounding());
    mpfr_add_d(result.mpfr_data(), result.mpfr_data(), value.x[3], REAL_REF::default_rounding());
    return result;
}
inline COMPLEX_REF cast2ref(const qd_complex &value) {
    return COMPLEX_REF(cast2ref(value.real()), cast2ref(value.imag()));
}
#elif defined ___MPLAPACK_BUILD_WITH_DD___
inline REAL_REF cast2ref(const dd_real &value) {
    REAL_REF result = REAL_REF::with_precision(REAL_REF::default_precision());
    mpfr_set_d(result.mpfr_data(), value.x[0], REAL_REF::default_rounding());
    mpfr_add_d(result.mpfr_data(), result.mpfr_data(), value.x[1], REAL_REF::default_rounding());
    return result;
}
inline COMPLEX_REF cast2ref(const dd_complex &value) {
    return COMPLEX_REF(cast2ref(value.real()), cast2ref(value.imag()));
}
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
inline REAL_REF cast2ref(double value) { return REAL_REF(value); }
inline COMPLEX_REF cast2ref(const complex<double> &value) {
    return COMPLEX_REF(REAL_REF(value.real()), REAL_REF(value.imag()));
}
#elif defined ___MPLAPACK_BUILD_WITH_BINARY80___
inline REAL_REF cast2ref(mplapack_binary80_t value) {
    return REAL_REF(mplapack::gmpfrxx_adapter::make_binary80_source(value));
}
inline COMPLEX_REF cast2ref(const complex<mplapack_binary80_t> &value) {
    return COMPLEX_REF(mplapack::gmpfrxx_adapter::make_binary80_complex_source(
        value.real(), value.imag()));
}
#elif defined ___MPLAPACK_BUILD_WITH_BINARY128___
inline REAL_REF cast2ref(mplapack_binary128_t value) {
    return REAL_REF(mplapack::gmpfrxx_adapter::make_binary128_source(value));
}
inline COMPLEX_REF cast2ref(const complex<mplapack_binary128_t> &value) {
    return COMPLEX_REF(mplapack::gmpfrxx_adapter::make_binary128_complex_source(
        value.real(), value.imag()));
}
#endif

// Generate in the active backend, then convert once to the reference type.
void set_random_number(REAL_REF &reference_value, REAL &backend_value);
void set_random_number(COMPLEX_REF &reference_value, COMPLEX &backend_value);
void set_random_number1to2(REAL_REF &reference_value, REAL &backend_value);
void set_random_number1to2(COMPLEX_REF &reference_value, COMPLEX &backend_value);
void set_random_number(INTEGER_REF &reference_value, INTEGER &backend_value);
void set_random_number1to2(INTEGER_REF &reference_value, INTEGER &backend_value);

REAL_REF infnorm(COMPLEX_REF *vec_ref, COMPLEX *vec, int len, int inc);
REAL_REF infnorm(REAL_REF *vec_ref, REAL *vec, int len, int inc);
INTEGER_REF infnorm(INTEGER_REF *vec_ref, INTEGER *vec, int len, int inc);
REAL_REF reldiff(REAL_REF *vec_ref, REAL *vec, int len, int inc);
REAL_REF infnorm_mat(REAL_REF *vec_ref, REAL *vec, int n, int m, int lda);

inline void set_random_vector(INTEGER_REF *vec_ref, INTEGER *vec, int len, int n) {
    if (len <= 0)
        return;
    for (int i = 0; i < len; i++) {
        vec_ref[i] = int(drand48() * (n - 1)) + 1;
        vec[i] = vec_ref[i];
    }
}

inline int veclen(int n, int inc) { return std::max(1, (1 + (abs(n) - 1) * abs(inc))); }

inline int vecplen(int n) {
    if (n <= 0)
        n = 1;
    return n * (n + 1) / 2;
}

inline int matlen(int lda, int n) { return std::max(1, abs(lda) * abs(n)); }

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
mpfr_class mpf_randomnumber(mpfr_class);
mpc_class mpc_randomnumber(mpc_class);
#endif
#if defined ___MPLAPACK_BUILD_WITH_GMP___
mpf_class mpf_randomnumber(mpf_class);
mpfc_class mpc_randomnumber(mpfc_class);
#endif
#if defined ___MPLAPACK_BUILD_WITH_QD___
qd_real mpf_randomnumber(qd_real);
qd_complex mpc_randomnumber(qd_complex);
#endif
#if defined ___MPLAPACK_BUILD_WITH_DD___
dd_real mpf_randomnumber(dd_real);
dd_complex mpc_randomnumber(dd_complex);
#endif
#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
double mpf_randomnumber(double);
complex<double> mpc_randomnumber(complex<double>);
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
mplapack_binary80_t mpf_randomnumber(mplapack_binary80_t dummy);
complex<mplapack_binary80_t> mpc_randomnumber(complex<mplapack_binary80_t> dummy);
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
mplapack_binary128_t mpf_randomnumber(mplapack_binary128_t dummy);
complex<mplapack_binary128_t> mpc_randomnumber(complex<mplapack_binary128_t> dummy);
#endif

template <class X_REF, class X> void set_random_vector(X_REF *vec_ref, X *vec, int len) {
    if (len <= 0)
        return;
    for (int i = 0; i < len; i++) {
        set_random_number(vec_ref[i], vec[i]);
    }
}

template <class X_REF, class X> void set_random_psdmat(X_REF *p_ref, X *p, int ldp, int n) {
    X_REF *tmpmat_ref = new X_REF[matlen(ldp, n)];
    X_REF mtmp_ref;
    X *tmpmat = new X[matlen(ldp, n)];
    X mtmp;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            set_random_number1to2(tmpmat_ref[i + j * ldp], tmpmat[i + j * ldp]);
        }
    }
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            mtmp_ref = 0.0;
            mtmp = 0.0;
            for (int k = 0; k < n; k++) {
                mtmp_ref = mtmp_ref + tmpmat_ref[i + k * ldp] * tmpmat_ref[j + k * ldp];
                mtmp = mtmp + tmpmat[i + k * ldp] * tmpmat[j + k * ldp];
            }
            p_ref[i + j * ldp] = mtmp_ref;
            p[i + j * ldp] = mtmp;
        }
    }
    delete[] tmpmat;
    delete[] tmpmat_ref;
}

template <class X_REF, class X> void set_random_symmmat_cond(X_REF *p_ref, X *p, int ldp, int n, int cond) {
    using cond_real = X;
    cond_real *tmpmat1 = new cond_real[matlen(ldp, n)];
    cond_real *tmpmat2 = new cond_real[matlen(ldp, n)];
    cond_real *tmpmat3 = new cond_real[matlen(ldp, n)];
    cond_real rtmp;

    for (int i = 0; i < matlen(ldp, n); i++) {
        p[i] = 0.0;
        p_ref[i] = 0.0;
    }
    for (int i = 0; i < n; i++) {
        rtmp = cond_real(cond) * cond_real(2.0 * (i + 1) - n) / cond_real(2.0 * n);
        tmpmat1[i + i * ldp] = pow(cond_real(10.0), rtmp);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmpmat2[i + j * ldp] = (mpf_randomnumber(rtmp) + cond_real(1.0)) / cond_real(2.0);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            rtmp = 0.0;
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    rtmp = rtmp + tmpmat2[i + k * ldp] * tmpmat1[k + l * ldp] * tmpmat2[j + l * ldp];
                }
            }
            tmpmat3[i + j * ldp] = rtmp;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            p[i + j * ldp] = tmpmat3[i + j * ldp];
            p_ref[i + j * ldp] = cast2ref(tmpmat3[i + j * ldp]);
        }
    }
    delete[] tmpmat1;
    delete[] tmpmat2;
    delete[] tmpmat3;
}

template <class X_REF, class X> void set_hilbertmat(X_REF *p_ref, X *p, int ldp, int n) {
    X_REF mtmp_ref, One_ref;
    X mtmp, One;
    One_ref = 1.0;
    One = 1.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mtmp_ref = (i + 1) + (j + 1) - 1;
            mtmp = (i + 1) + (j + 1) - 1;
            p_ref[i + j * ldp] = One_ref / mtmp_ref;
            p[i + j * ldp] = One / mtmp;
        }
    }
}

#include <mplapack_utils.h>

#endif

/*
 * Copyright (c) 2021
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

#ifndef _MPLAPACK_DEBUG_H_
#define _MPLAPACK_DEBUG_H_

#include <cmath>
using namespace std;

#define GMP_FORMAT            "%+36.32Fe"
#define GMP_SHORT_FORMAT      "%+20.16Fe"
#define MPFR_FORMAT           "%+68.64Re"
#define MPFR_SHORT_FORMAT     "%+20.16Re"
#define DOUBLE_FORMAT         "%+20.16e"
#define DOUBLE_SHORT_FORMAT   "%+20.16e"
#define FLOAT64X_FORMAT       "%+24.20Le"
#define FLOAT64X_SHORT_FORMAT "%+20.16Le"

#define FLOAT128_FORMAT       "%+24.20Q"
#define FLOAT128_SHORT_FORMAT "%+20.16Q"

#define QD_PRECISION         64
#define QD_PRECISION_SHORT   16
#define DD_PRECISION         32
#define DD_PRECISION_SHORT   16
#define BUFLEN 1024

void printnum(double rtmp);
void printnum(complex<double> ctmp);
void printnum(long double rtmp);
void printnum(complex<long double> ctmp);
void printnum(__complex__ double ctmp);

inline void printnum_short(INTEGER itmp) {
    printf("%d ", (int)itmp);
}

void printnum_short(double rtmp);
void printnum_short(complex<double> ctmp);
void printnum_short(__complex__ double ctmp);

void sprintnum(char *buf, double rtmp);
void sprintnum(char *buf, complex<double> rtmp);
void sprintnum(char *buf, __complex__ double rtmp);

void sprintnum_short(char *buf, double rtmp);
void sprintnum_short(char *buf, complex<double> rtmp);
void sprintnum_short(char *buf, __complex__ double rtmp);

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
inline void printnum(mpreal rtmp) {
    mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void printnum_short(mpreal rtmp) {
    mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void printnum(mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_FORMAT MPFR_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
inline void printnum_short(mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
inline void sprintnum(char *buf, mpreal rtmp) {
    mpfr_snprintf(buf, BUFLEN, MPFR_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum_short(char *buf, mpreal rtmp) {
    mpfr_snprintf(buf, BUFLEN, MPFR_SHORT_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum(char *buf, mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, BUFLEN, MPFR_FORMAT MPFR_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
inline void sprintnum_short(char *buf, mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, BUFLEN, MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void printnum(mpf_class rtmp);
void printnum(mpc_class ctmp);
void printnum_short(mpf_class rtmp);
void printnum_short(mpc_class ctmp);
void sprintnum(char *buf, mpf_class rtmp);
void sprintnum(char *buf, mpc_class rtmp);
void sprintnum_short(char *buf, mpf_class rtmp);
void sprintnum_short(char *buf, mpc_class rtmp);
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
void printnum(qd_real rtmp);
void printnum(qd_complex rtmp);
void printnum_short(qd_real rtmp);
void printnum_short(qd_complex rtmp);
void sprintnum(char *buf, qd_real rtmp);
void sprintnum(char *buf, qd_complex rtmp);
void sprintnum_short(char *buf, qd_real rtmp);
void sprintnum_short(char *buf, qd_complex rtmp);
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
void printnum(dd_real rtmp);
void printnum(dd_complex rtmp);
void printnum_short(dd_real rtmp);
void printnum_short(dd_complex rtmp);
void sprintnum(char *buf, dd_real rtmp);
void sprintnum(char *buf, dd_complex rtmp);
void sprintnum_short(char *buf, dd_real rtmp);
void sprintnum_short(char *buf, dd_complex rtmp);
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
void printnum(_Float64x rtmp);
void printnum(complex<_Float64x> rtmp);
void printnum_short(_Float64x rtmp);
void printnum_short(complex<_Float64x> rtmp);
void sprintnum(char *buf, _Float64x rtmp);
void sprintnum(char *buf, complex<_Float64x> rtmp);
void sprintnum_short(char *buf, _Float64x rtmp);
void sprintnum_short(char *buf, complex<_Float64x> rtmp);
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
void printnum(_Float128 rtmp);
void printnum(complex<_Float128> rtmp);
void printnum_short(_Float128 rtmp);
void printnum_short(complex<_Float128> rtmp);
void sprintnum(char *buf, _Float128 rtmp);
void sprintnum(char *buf, complex<_Float128> rtmp);
void sprintnum_short(char *buf, _Float128 rtmp);
void sprintnum_short(char *buf, complex<_Float128> rtmp);
#endif

template <class X> void printmat(int N, int M, X *A, int LDA) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("[ ");
        for (int j = 0; j < M; j++) {
            tmp = A[i + j * LDA];
            printnum_short(tmp);
            if (j < M - 1)
                printf(", ");
        }
        if (i < N - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}

template <class X> void printmatU(int N, X *A, int LDA) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("[ ");
        for (int j = 0; j < N; j++) {
	  if (i<=j)
            tmp = A[i + j * LDA];
	  else
            tmp = A[j + i * LDA];
            printnum_short(tmp);
            if (j < N - 1)
                printf(", ");
        }
        if (i < N - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}

template <class X> void printvec(X *A, int len) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = A[i];
        printnum_short(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

#if defined __MPLAPACK_XLAENV__
#define _MPLAPACK_XLAENV_EXTERN_
#else
#define _MPLAPACK_XLAENV_EXTERN_ extern
#endif

_MPLAPACK_XLAENV_EXTERN_ INTEGER iparms[100];

#if defined __MPLAPACK_MXERBLA__
#define _MPLAPACK_MXERBLA_EXTERN_
#else
#define _MPLAPACK_MXERBLA_EXTERN_ extern
#endif
#define srnamt_len 1024
_MPLAPACK_MXERBLA_EXTERN_ int infot;
_MPLAPACK_MXERBLA_EXTERN_ INTEGER nout;
_MPLAPACK_MXERBLA_EXTERN_ bool ok;
_MPLAPACK_MXERBLA_EXTERN_ bool lerr;
_MPLAPACK_MXERBLA_EXTERN_ char srnamt[srnamt_len];

#endif

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

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }
inline void printnum(std::complex<double> ctmp) { printf(DOUBLE_FORMAT DOUBLE_FORMAT "i", ctmp.real(), ctmp.imag()); }
inline void printnum(__complex__ double ctmp) { printf(DOUBLE_FORMAT DOUBLE_FORMAT "i", __real__ ctmp, __imag__ ctmp); }

inline void printnum_short(double rtmp) { printf(DOUBLE_SHORT_FORMAT, rtmp); }
inline void printnum_short(std::complex<double> ctmp) { printf(DOUBLE_SHORT_FORMAT DOUBLE_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }
inline void printnum_short(__complex__ double ctmp) { printf(DOUBLE_SHORT_FORMAT DOUBLE_SHORT_FORMAT "i", __real__ ctmp, __imag__ ctmp); }
inline void printnum_short(mplapackint itmp) { printf("%d ", (int)itmp); }

inline void sprintnum(char *buf, double rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, DOUBLE_FORMAT, rtmp); }
inline void sprintnum(char *buf, std::complex<double> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, DOUBLE_FORMAT DOUBLE_FORMAT "i", ctmp.real(), ctmp.imag()); }
inline void sprintnum(char *buf, __complex__ double ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, DOUBLE_FORMAT DOUBLE_FORMAT "i", __real__ ctmp, __imag__ ctmp); }

inline void sprintnum_short(char *buf, double rtmp) { snprintf(buf, __MPLAPACK_BUFLEN__, DOUBLE_SHORT_FORMAT, rtmp); }
inline void sprintnum_short(char *buf, std::complex<double> ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, DOUBLE_SHORT_FORMAT DOUBLE_SHORT_FORMAT "i", ctmp.real(), ctmp.imag()); }
inline void sprintnum_short(char *buf, __complex__ double ctmp) { snprintf(buf, __MPLAPACK_BUFLEN__, DOUBLE_SHORT_FORMAT DOUBLE_SHORT_FORMAT "i", __real__ ctmp, __imag__ ctmp); }

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
            if (i <= j)
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

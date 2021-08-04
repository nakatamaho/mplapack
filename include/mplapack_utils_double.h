/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_utils_double.h,v 1.11 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MUTILS_DOUBLE_H_
#define _MUTILS_DOUBLE_H_

#if !(_XOPEN_SOURCE >= 600) && !(_ISOC99_SOURCE)
inline double log2(double x) { return log(x) / log(2.0); }
#endif

#include <complex>

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

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

inline double pow2(double a) {
    double mtmp;
    mtmp = a * a;
    return mtmp;
}

// implementation of sign transfer function.
inline double sign(double a, double b) {
    double mtmp;
    mtmp = std::abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline double castREAL_double(mplapackint n) {
    double ret = n;
    return ret;
}

inline mplapackint castINTEGER_double(double a) {
    mplapackint i = a;
    return i;
}

inline long nint(double a) {
    long i;
    double tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (long)tmp;
    return i;
}

inline double cast2double(double a) { return a; }

inline std::complex<double> exp(std::complex<double> x) {
    double ex;
    double c;
    double s;
    std::complex<double> ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

inline double pi(double dummy) { return M_PI; }

#endif

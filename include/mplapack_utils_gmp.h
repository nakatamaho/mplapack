/*
 * Copyright (c) 2008-2021
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

#include "mpc_class.h"

#define GMP_FORMAT "%+36.32Fe"
#define GMP_SHORT_FORMAT "%+20.16Fe"
#define BUFLEN 1024

inline void printnum(mpf_class rtmp) {
    gmp_printf(GMP_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void printnum_short(mpf_class rtmp) {
    gmp_printf(GMP_SHORT_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void printnum(mpc_class ctmp) {
    gmp_printf(GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void printnum_short(mpc_class ctmp) {
    gmp_printf(GMP_SHORT_FORMAT GMP_SHORT_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void sprintnum(char *buf, mpf_class rtmp) {
    gmp_snprintf(buf, BUFLEN, GMP_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void sprintnum_short(char *buf, mpf_class rtmp) {
    gmp_snprintf(buf, BUFLEN, GMP_SHORT_FORMAT, rtmp.get_mpf_t());
    return;
}

inline void sprintnum(char *buf, mpc_class ctmp) {
    gmp_snprintf(buf, BUFLEN, GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline void sprintnum_short(char *buf, mpc_class ctmp) {
    gmp_snprintf(buf, BUFLEN, GMP_SHORT_FORMAT GMP_SHORT_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

inline mpf_class sign(mpf_class a, mpf_class b) {
    mpf_class mtmp;
    mpf_abs(mtmp.get_mpf_t(), a.get_mpf_t());
    if (b != 0.0) {
        mtmp = mpf_sgn(b.get_mpf_t()) * mtmp;
    }
    return mtmp;
}

inline mpf_class atan2(mpf_class a, mpf_class b) {
    double dtemp1, dtemp2;
    mpf_class mtemp3;
    dtemp1 = a.get_d();
    dtemp2 = b.get_d();
    mtemp3 = mpf_class(atan2(a, b));
    return mtemp3;
}

inline mpf_class castREAL_gmp(mplapackint n) {
    mpf_class a(n);
    return a;
}

inline mplapackint castINTEGER_gmp(mpf_class a) {
    mplapackint i;
    i = mpf_get_si(a.get_mpf_t());
    return i;
}

inline mplapackint nint(mpf_class a) {
    mplapackint i;
    mpf_class tmp;
    a = a + 0.5;
    mpf_floor(tmp.get_mpf_t(), a.get_mpf_t());
    i = mpf_get_si(tmp.get_mpf_t());
    return i;
}

inline double cast2double(mpf_class a) { return a.get_d(); }

inline mpc_class sin(mpc_class a) {
    double dtemp1, dtemp2;
    mpf_class mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    dtemp1 = mtemp1.get_d();
    dtemp2 = mtemp2.get_d();
    mtemp1 = sin(dtemp1) * cosh(dtemp2);
    mtemp2 = cos(dtemp1) * sinh(dtemp2);
    mpc_class b = mpc_class(mtemp1, mtemp2);
    return b;
}

inline mpf_class log2(mpf_class x) {
    double d;
    double ln2_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln2_app = (double)exp + log10(d) / log10(2);
    return ln2_app;
}

inline mpf_class log(mpf_class x) {
    double d;
    double ln_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln_app = (double)exp * log(2.0) + log(d);
    return ln_app;
}

inline mpf_class log10(mpf_class x) {
    double d;
    double ln10_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln10_app = (double)exp * log10(2.0) + log10(d);
    return ln10_app;
}

inline mpf_class pow(mpf_class x, mpf_class y) {
    mpf_class mtemp1, mtemp2;
    mtemp1 = y * log(x);
    mtemp2 = exp(mtemp1.get_d());
    return mtemp2;
}

inline mpf_class cos(mpf_class x) {
    mpf_class mtemp1;
    mtemp1 = cos(x.get_d());
    return mtemp1;
}

inline mpf_class sin(mpf_class x) {
    mpf_class mtemp1;
    mtemp1 = sin(x.get_d());
    return mtemp1;
}

inline mpf_class exp(mpf_class x) {
    mpf_class mtemp1;
    mtemp1 = exp(x.get_d());
    return mtemp1;
}

inline mpf_class pi(mpf_class dummy) {
    mpf_class mtemp1;
    mtemp1 = M_PI; // returns pi in double(!)
    return mtemp1;
}

inline mpc_class exp(mpc_class x) {
    mpf_class ex;
    mpf_class c;
    mpf_class s;
    mpc_class ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

#endif

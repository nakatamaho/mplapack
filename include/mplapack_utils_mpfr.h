/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_utils_mpfr.h,v 1.9 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MUTILS_MPFR_H_
#define _MUTILS_MPFR_H_

#include "mpcomplex.h"
#include "mpreal.h"

using namespace mpfr;

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

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
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum_short(char *buf, mpreal rtmp) {
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_SHORT_FORMAT, mpfr_ptr(rtmp));
    return;
}
inline void sprintnum(char *buf, mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_FORMAT MPFR_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}
inline void sprintnum_short(char *buf, mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, __MPLAPACK_BUFLEN__, MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

inline mpreal pow2(mpreal a) {
    mpreal mtmp = a * a;
    return mtmp;
}

// implementation of sign transfer function.
inline mpreal sign(mpreal a, mpreal b) {
    mpreal mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline mplapackint nint(mpreal a) {
    mplapackint i;
    mpreal tmp;
    a = a + 0.5;
    tmp = floorl(a);
    i = tmp; // cast to long
    return i;
}

inline mplapackint castINTEGER_mpfr(mpreal a) {
    mplapackint i;
    i = a;
    return i;
}

inline mpreal castREAL_mpfr(mplapackint a) {
    mpreal i = a;
    return i;
}

inline mpreal pi(mpreal dummy) {
    mpreal _PI;
    _PI = const_pi(mpreal::get_default_prec());
    return _PI;
}

#endif

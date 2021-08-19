/*
 * Copyright (c) 2008-2021
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

#ifndef _MUTILS_DD_H_
#define _MUTILS_DD_H_

#if defined ___MPLAPACK_INTERNAL___
#define DD_PRECISION 32
#define DD_PRECISION_SHORT 16

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

#include <cstring>

inline void printnum(dd_real rtmp) {
    std::cout.precision(DD_PRECISION);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

inline void printnum_short(dd_real rtmp) {
    std::cout.precision(DD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}
inline void printnum(dd_complex rtmp) {
    std::cout.precision(DD_PRECISION);
    if (rtmp.real() >= 0.0) {
        std::cout << "+" << rtmp.real();
    } else {
        std::cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
        std::cout << "+" << rtmp.imag() << "i";
    } else {
        std::cout << rtmp.imag() << "i";
    }
    return;
}
inline void printnum_short(dd_complex rtmp) {
    std::cout.precision(DD_PRECISION);
    if (rtmp.real() >= 0.0) {
        std::cout << "+" << rtmp.real();
    } else {
        std::cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
        std::cout << "+" << rtmp.imag() << "i";
    } else {
        std::cout << rtmp.imag() << "i";
    }
    return;
}
inline void sprintnum(char *buf, dd_real rtmp) {
    rtmp.write(buf, __MPLAPACK_BUFLEN__, DD_PRECISION);
    return;
}
inline void sprintnum_short(char *buf, dd_real rtmp) {
    rtmp.write(buf, __MPLAPACK_BUFLEN__, DD_PRECISION_SHORT);
    return;
}
inline void sprintnum(char *buf, dd_complex rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    rtmp.real().write(buf1, __MPLAPACK_BUFLEN__, DD_PRECISION);
    rtmp.real().write(buf2, __MPLAPACK_BUFLEN__, DD_PRECISION);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
inline void sprintnum_short(char *buf, dd_complex rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    rtmp.real().write(buf1, __MPLAPACK_BUFLEN__, DD_PRECISION_SHORT);
    rtmp.real().write(buf2, __MPLAPACK_BUFLEN__, DD_PRECISION_SHORT);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
#endif

inline dd_real pow2(dd_real a) {
    dd_real mtmp = a * a;
    return mtmp;
}

// implementation of sign transfer function.
inline dd_real sign(dd_real a, dd_real b) {
    dd_real mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline dd_real castREAL_dd(mplapackint n) {
    dd_real ret;
    ret.x[0] = (static_cast<double>(n));
    ret.x[1] = 0.0;
    return ret;
}
inline mplapackint castINTEGER_dd(dd_real a) {
    mplapackint i = a.x[0];
    return i;
}

inline long __dd_nint(dd_real a) {
    long i;
    dd_real tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (int)tmp.x[0];
    return i;
}

inline double cast2double(dd_real a) { return a.x[0]; }

inline dd_complex sin(dd_complex a) {
    dd_real mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    dd_complex b = dd_complex(sin(mtemp1) * cosh(mtemp2), cos(mtemp1) * sinh(mtemp2));
    return b;
}

inline dd_complex cos(dd_complex a) {
    dd_real mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    dd_complex b = dd_complex(cos(mtemp1) * cosh(mtemp2), -sin(mtemp1) * sinh(mtemp2));
    return b;
}

inline dd_real log2(dd_real x) { return log10(x) / (dd_real::_log2 / dd_real::_log10); }

inline dd_complex exp(dd_complex x) {
    dd_real ex;
    dd_real c;
    dd_real s;
    dd_complex ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

inline dd_real pi(dd_real dummy) { return dd_real::_pi; }

#endif

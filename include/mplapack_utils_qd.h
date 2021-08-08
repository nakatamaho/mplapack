/*
 * Copyright (c) 2008-2021
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

#ifndef _MUTILS_QD_H_
#define _MUTILS_QD_H_

#if defined ___MPLAPACK_INTERNAL___
#define QD_PRECISION 64
#define QD_PRECISION_SHORT 16

#include <cstring>

#if !defined __MPLAPACK_BUFLEN__
#define __MPLAPACK_BUFLEN__ 1024
#endif

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

inline void printnum_short(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

inline void printnum(qd_complex rtmp) {
    std::cout.precision(QD_PRECISION);
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

inline void printnum_short(qd_complex rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
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

inline void sprintnum(char *buf, qd_real rtmp) {
    rtmp.write(buf, __MPLAPACK_BUFLEN__, QD_PRECISION);
    return;
}
inline void sprintnum_short(char *buf, qd_real rtmp) {
    rtmp.write(buf, __MPLAPACK_BUFLEN__, QD_PRECISION_SHORT);
    return;
}
inline void sprintnum(char *buf, qd_complex rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    rtmp.real().write(buf1, __MPLAPACK_BUFLEN__, QD_PRECISION);
    rtmp.real().write(buf2, __MPLAPACK_BUFLEN__, QD_PRECISION);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
inline void sprintnum_short(char *buf, qd_complex rtmp) {
    char buf1[__MPLAPACK_BUFLEN__], buf2[__MPLAPACK_BUFLEN__];
    rtmp.real().write(buf1, __MPLAPACK_BUFLEN__, QD_PRECISION_SHORT);
    rtmp.real().write(buf2, __MPLAPACK_BUFLEN__, QD_PRECISION_SHORT);
    strcat(buf, buf1);
    strcat(buf, buf2);
    strcat(buf, "i");
}
#endif

inline qd_real pow2(qd_real a) {
    qd_real mtmp = a * a;
    return mtmp;
}

// implementation of sign transfer function.
inline qd_real sign(qd_real a, qd_real b) {
    qd_real mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
        mtmp = -mtmp;
    }
    return mtmp;
}

inline qd_real castREAL_qd(mplapackint n) {
    qd_real ret;
    ret.x[0] = (static_cast<double>(n));
    ret.x[1] = 0.0;
    return ret;
}
inline mplapackint castINTEGER_qd(qd_real a) {
    mplapackint i = a.x[0];
    return i;
}

inline long __qd_nint(qd_real a) {
    long i;
    qd_real tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (int)tmp.x[0];
    return i;
}

inline double cast2double(qd_real a) { return a.x[0]; }

inline qd_complex sin(qd_complex a) {
    qd_real mtemp1, mtemp2;
    mtemp1 = a.real();
    mtemp2 = a.imag();
    qd_complex b = qd_complex(sin(mtemp1) * cosh(mtemp2), cos(mtemp1) * sinh(mtemp2));
    return b;
}

inline qd_real log2(qd_real x) { return log10(x) / (qd_real::_log2 / qd_real::_log10); }

inline qd_complex exp(qd_complex x) {
    qd_real ex;
    qd_real c;
    qd_real s;
    qd_complex ans;
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
    ans.real(ex * c);
    ans.imag(ex * s);
    return ans;
}

inline qd_real pi(qd_real dummy) { return qd_real::_pi; }
#endif

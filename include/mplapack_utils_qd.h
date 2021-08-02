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

qd_real log2(qd_real x);
qd_complex exp(qd_complex x);
qd_real pi(qd_real dummy);

qd_real sign(qd_real a, qd_real b);
qd_complex Real2Complex(qd_real a, qd_real b);
qd_complex sin(qd_complex a);

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

inline qd_complex Real2Complex(qd_real a, qd_real b) {
    qd_complex ret;
    ret.real(a);
    ret.imag(b);
    return ret;
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

#endif

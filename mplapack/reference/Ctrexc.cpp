/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
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

#include <mpblas.h>
#include <mplapack.h>

void Ctrexc(const char *compq, INTEGER const n, COMPLEX *t, INTEGER const ldt, COMPLEX *q, INTEGER const ldq, INTEGER const ifst, INTEGER const ilst, INTEGER &info) {
    //
    //     Decode and test the input parameters.
    //
    info = 0;
    bool wantq = Mlsame(compq, "V");
    if (!Mlsame(compq, "N") && !wantq) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (ldt < max((INTEGER)1, n)) {
        info = -4;
    } else if (ldq < 1 || (wantq && ldq < max((INTEGER)1, n))) {
        info = -6;
    } else if ((ifst < 1 || ifst > n) && (n > 0)) {
        info = -7;
    } else if ((ilst < 1 || ilst > n) && (n > 0)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Ctrexc", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1 || ifst == ilst) {
        return;
    }
    //
    INTEGER m1 = 0;
    INTEGER m2 = 0;
    INTEGER m3 = 0;
    if (ifst < ilst) {
        //
        //        Move the IFST-th diagonal element forward down the diagonal.
        //
        m1 = 0;
        m2 = -1;
        m3 = 1;
    } else {
        //
        //        Move the IFST-th diagonal element backward up the diagonal.
        //
        m1 = -1;
        m2 = 0;
        m3 = -1;
    }
    //
    INTEGER k = 0;
    COMPLEX t11 = 0.0;
    COMPLEX t22 = 0.0;
    REAL cs = 0.0;
    COMPLEX sn = 0.0;
    COMPLEX temp = 0.0;
    for (k = ifst + m1; m3 >=0 ? k <= ilst + m2 : k >= ilst + m2 ; k = k + m3) {
        //
        //        Interchange the k-th and (k+1)-th diagonal elements.
        //
        t11 = t[(k - 1) + (k - 1) * ldt];
        t22 = t[((k + 1) - 1) + ((k + 1) - 1) * ldt];
        //
        //        Determine the transformation to perform the interchange.
        //
        Clartg(t[(k - 1) + ((k + 1) - 1) * ldt], t22 - t11, cs, sn, temp);
        //
        //        Apply transformation to the matrix T.
        //
        if (k + 2 <= n) {
            Crot(n - k - 1, &t[(k - 1) + ((k + 2) - 1) * ldt], ldt, &t[((k + 1) - 1) + ((k + 2) - 1) * ldt], ldt, cs, sn);
        }
        Crot(k - 1, &t[(k - 1) * ldt], 1, &t[((k + 1) - 1) * ldt], 1, cs, conj(sn));
        //
        t[(k - 1) + (k - 1) * ldt] = t22;
        t[((k + 1) - 1) + ((k + 1) - 1) * ldt] = t11;
        //
        if (wantq) {
            //
            //           Accumulate transformation in the matrix Q.
            //
            Crot(n, &q[(k - 1) * ldq], 1, &q[((k + 1) - 1) * ldq], 1, cs, conj(sn));
        }
        //
    }
    //
    //     End of Ctrexc
    //
}

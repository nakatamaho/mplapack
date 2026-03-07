/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZLARF1L.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Clarf1l(const char *side, INTEGER const m, INTEGER const n, COMPLEX *v, INTEGER const incv, COMPLEX const tau, COMPLEX *c, INTEGER const ldc, COMPLEX *work) {
    //
    bool applyleft = Mlsame(side, "L");
    INTEGER firstv = 1;
    INTEGER lastc = 0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    INTEGER lastv = 0;
    INTEGER i = 0;
    if (tau != zero) {
        // Set up variables for scanning V.  LASTV begins pointing to the end
        // of V up to V(1).
        if (applyleft) {
            lastv = m;
        } else {
            lastv = n;
        }
        i = 1;
        // Look for the last non-zero row in V.
        while (lastv > firstv && v[i - 1] == zero) {
            firstv++;
            i += incv;
        }
        if (applyleft) {
            // Scan for the last non-zero column in C(1:lastv,:).
            lastc = iMlazlc(lastv, n, c, ldc);
        } else {
            // Scan for the last non-zero row in C(:,1:lastv).
            lastc = iMlazlr(m, lastv, c, ldc);
        }
    }
    if (lastc == 0) {
        return;
    }
    const COMPLEX one = COMPLEX(1.0, 0.0);
    INTEGER j = 0;
    if (applyleft) {
        //
        // Form  H * C
        //
        if (lastv == firstv) {
            //
            // C(lastv,1:lastc) := ( 1 - tau ) * C(lastv,1:lastc)
            //
            Cscal(lastc, one - tau, &c[(lastv - 1)], ldc);
        } else {
            //
            // w(1:lastc,1) := C(firstv:lastv-1,1:lastc)**T * v(firstv:lastv-1,1)
            //
            Cgemv("Conjugate transpose", lastv - firstv, lastc, one, &c[(firstv - 1)], ldc, &v[i - 1], incv, zero, work, 1);
            //
            // w(1:lastc,1) += C(lastv,1:lastc)**H * v(lastv,1)
            //
            for (j = 1; j <= lastc; j = j + 1) {
                work[j - 1] += conj(c[(lastv - 1) + (j - 1) * ldc]);
            }
            //
            // C(lastv,1:lastc) += - tau * v(lastv,1) * w(1:lastc,1)**H
            //
            for (j = 1; j <= lastc; j = j + 1) {
                c[(lastv - 1) + (j - 1) * ldc] = c[(lastv - 1) + (j - 1) * ldc] - tau * conj(work[j - 1]);
            }
            //
            // C(firstv:lastv-1,1:lastc) += - tau * v(firstv:lastv-1,1) * w(1:lastc,1)**H
            //
            Cgerc(lastv - firstv, lastc, -tau, &v[i - 1], incv, work, 1, &c[(firstv - 1)], ldc);
        }
    } else {
        //
        // Form  C * H
        //
        if (lastv == firstv) {
            //
            // C(1:lastc,lastv) := ( 1 - tau ) * C(1:lastc,lastv)
            //
            Cscal(lastc, one - tau, &c[(lastv - 1) * ldc], 1);
        } else {
            //
            // w(1:lastc,1) := C(1:lastc,firstv:lastv-1) * v(firstv:lastv-1,1)
            //
            Cgemv("No transpose", lastc, lastv - firstv, one, &c[(firstv - 1) * ldc], ldc, &v[i - 1], incv, zero, work, 1);
            //
            // w(1:lastc,1) += C(1:lastc,lastv) * v(lastv,1)
            //
            Caxpy(lastc, one, &c[(lastv - 1) * ldc], 1, work, 1);
            //
            // C(1:lastc,lastv) += - tau * v(lastv,1) * w(1:lastc,1)
            //
            Caxpy(lastc, -tau, work, 1, &c[(lastv - 1) * ldc], 1);
            //
            // C(1:lastc,firstv:lastv-1) += - tau * w(1:lastc,1) * v(firstv:lastv-1)**H
            //
            Cgerc(lastc, lastv - firstv, -tau, work, 1, &v[i - 1], incv, &c[(firstv - 1) * ldc], ldc);
        }
    }
    //
    // End of Clarf1l
    //
}

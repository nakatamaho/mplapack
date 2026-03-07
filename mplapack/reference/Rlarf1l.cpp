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

// Derived from LAPACK routine DLARF1L.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Rlarf1l(const char *side, INTEGER const m, INTEGER const n, REAL *v, INTEGER const incv, REAL const tau, REAL *c, INTEGER const ldc, REAL *work) {
    //
    bool applyleft = Mlsame(side, "L");
    INTEGER firstv = 1;
    INTEGER lastc = 0;
    const REAL zero = 0.0;
    INTEGER lastv = 0;
    INTEGER i = 0;
    if (tau != zero) {
        // Set up variables for scanning V.  LASTV begins pointing to the end
        // of V.
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
            lastc = iMladlc(lastv, n, c, ldc);
        } else {
            // Scan for the last non-zero row in C(:,1:lastv).
            lastc = iMladlr(m, lastv, c, ldc);
        }
    }
    if (lastc == 0) {
        return;
    }
    const REAL one = 1.0;
    if (applyleft) {
        //
        // Form  H * C
        //
        if (lastv > 0) {
            //
            // C := HC = (1-\tau)C.
            if (lastv == firstv) {
                Rscal(lastc, one - tau, &c[(firstv - 1)], ldc);
            } else {
                //
                // w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
                //
                Rgemv("Transpose", lastv - firstv, lastc, one, &c[(firstv - 1)], ldc, &v[i - 1], incv, zero, work, 1);
                //
                Raxpy(lastc, one, &c[(lastv - 1)], ldc, work, 1);
                //
                // C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T
                //
                // = C(...) - tau * w(1:lastc,1)**T
                Raxpy(lastc, -tau, work, 1, &c[(lastv - 1)], ldc);
                //
                Rger(lastv - firstv, lastc, -tau, &v[i - 1], incv, work, 1, &c[(firstv - 1)], ldc);
            }
        }
    } else {
        //
        // Form  C * H
        //
        if (lastv > 0) {
            //
            // C := CH = C(1-\tau).
            if (lastv == firstv) {
                Rscal(lastc, one - tau, c, 1);
            } else {
                //
                // w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
                //
                // w(1:lastc,1) := C(1:lastc,1:lastv-1) * v(1:lastv-1,1)
                Rgemv("No transpose", lastc, lastv - firstv, one, &c[(firstv - 1) * ldc], ldc, &v[i - 1], incv, zero, work, 1);
                //
                Raxpy(lastc, one, &c[(lastv - 1) * ldc], 1, work, 1);
                //
                // C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T
                //
                // = C(...) - tau * w(1:lastc,1)
                Raxpy(lastc, -tau, work, 1, &c[(lastv - 1) * ldc], 1);
                //
                Rger(lastc, lastv - firstv, -tau, work, 1, &v[i - 1], incv, &c[(firstv - 1) * ldc], ldc);
            }
        }
    }
    //
    // End of Rlarf1l
    //
}

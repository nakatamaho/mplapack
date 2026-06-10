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

// Derived from LAPACK routine DLAED3.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Rlaed3(INTEGER const k, INTEGER const n, INTEGER const n1, REAL *d, REAL *q, INTEGER const ldq, REAL const rho, REAL *dlambda, REAL *q2, INTEGER *indx, INTEGER *ctot, REAL *w, REAL *s, INTEGER &info) {
    INTEGER j = 0;
    INTEGER ii = 0;
    INTEGER i = 0;
    REAL temp = 0.0;
    INTEGER n2 = 0;
    INTEGER n12 = 0;
    INTEGER n23 = 0;
    INTEGER iq2 = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    //
    // Test the input parameters.
    //
    info = 0;
    //
    if (k < 0) {
        info = -1;
    } else if (n < k) {
        info = -2;
    } else if (ldq < max((INTEGER)1, n)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rlaed3", -info);
        return;
    }
    //
    // Quick return if possible
    //
    if (k == 0) {
        return;
    }
    //
    for (j = 1; j <= k; j = j + 1) {
        Rlaed4(k, j, dlambda, w, &q[(j - 1) * ldq], rho, d[j - 1], info);
        //
        // If the zero finder fails, the computation is terminated.
        //
        if (info != 0) {
            goto statement_120;
        }
    }
    //
    if (k == 1) {
        goto statement_110;
    }
    if (k == 2) {
        for (j = 1; j <= k; j = j + 1) {
            w[1 - 1] = q[(j - 1) * ldq];
            w[2 - 1] = q[(2 - 1) + (j - 1) * ldq];
            ii = indx[1 - 1];
            q[(j - 1) * ldq] = w[ii - 1];
            ii = indx[2 - 1];
            q[(2 - 1) + (j - 1) * ldq] = w[ii - 1];
        }
        goto statement_110;
    }
    //
    // Compute updated W.
    //
    Rcopy(k, w, 1, s, 1);
    //
    // Initialize W(I) = Q(I,I)
    //
    Rcopy(k, q, ldq + 1, w, 1);
    for (j = 1; j <= k; j = j + 1) {
        for (i = 1; i <= j - 1; i = i + 1) {
            w[i - 1] = w[i - 1] * (q[(i - 1) + (j - 1) * ldq] / (dlambda[i - 1] - dlambda[j - 1]));
        }
        for (i = j + 1; i <= k; i = i + 1) {
            w[i - 1] = w[i - 1] * (q[(i - 1) + (j - 1) * ldq] / (dlambda[i - 1] - dlambda[j - 1]));
        }
    }
    for (i = 1; i <= k; i = i + 1) {
        w[i - 1] = sign(sqrt(-w[i - 1]), s[i - 1]);
    }
    //
    // Compute eigenvectors of the modified rank-1 modification.
    //
    for (j = 1; j <= k; j = j + 1) {
        for (i = 1; i <= k; i = i + 1) {
            s[i - 1] = w[i - 1] / q[(i - 1) + (j - 1) * ldq];
        }
        temp = Rnrm2(k, s, 1);
        for (i = 1; i <= k; i = i + 1) {
            ii = indx[i - 1];
            q[(i - 1) + (j - 1) * ldq] = s[ii - 1] / temp;
        }
    }
//
// Compute the updated eigenvectors.
//
statement_110:
    //
    n2 = n - n1;
    n12 = ctot[1 - 1] + ctot[2 - 1];
    n23 = ctot[2 - 1] + ctot[3 - 1];
    //
    Rlacpy("A", n23, k, &q[((ctot[1 - 1] + 1) - 1)], ldq, s, n23);
    iq2 = n1 * n12 + 1;
    if (n23 != 0) {
        Rgemm("N", "N", n2, k, n23, one, &q2[iq2 - 1], n2, s, n23, zero, &q[((n1 + 1) - 1)], ldq);
    } else {
        Rlaset("A", n2, k, zero, zero, &q[((n1 + 1) - 1)], ldq);
    }
    //
    Rlacpy("A", n12, k, q, ldq, s, n12);
    if (n12 != 0) {
        Rgemm("N", "N", n1, k, n12, one, q2, n1, s, n12, zero, q, ldq);
    } else {
        Rlaset("A", n1, k, zero, zero, &q[0], ldq);
    }
//
statement_120:;
    //
    // End of Rlaed3
    //
}

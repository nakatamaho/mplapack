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

inline REAL Rlamc3(REAL a, REAL b) {
    REAL c = a + b;
    return c;
}

void Rlasd8(INTEGER const icompq, INTEGER const k, REAL *d, REAL *z, REAL *vf, REAL *vl, REAL *difl, REAL *difr, INTEGER const lddifr, REAL *dsigma, REAL *work, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    info = 0;
    //
    if ((icompq < 0) || (icompq > 1)) {
        info = -1;
    } else if (k < 1) {
        info = -2;
    } else if (lddifr < k) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Rlasd8", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    const REAL one = 1.0;
    if (k == 1) {
        d[1 - 1] = abs(z[1 - 1]);
        difl[1 - 1] = d[1 - 1];
        if (icompq == 1) {
            difl[2 - 1] = one;
            difr[(2 - 1) * lddifr] = one;
        }
        return;
    }
    //
    //     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
    //     be computed with high relative accuracy (barring over/underflow).
    //     This is a problem on machines without a guard digit in
    //     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
    //     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
    //     which on any of these machines zeros out the bottommost
    //     bit of DSIGMA(I) if it is 1; this makes the subsequent
    //     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
    //     occurs. On binary machines with a guard digit (almost all
    //     machines) it does not change DSIGMA(I) at all. On hexadecimal
    //     and decimal machines with a guard digit, it slightly
    //     changes the bottommost bits of DSIGMA(I). It does not account
    //     for hexadecimal or decimal machines without guard digits
    //     (we know of none). We use a subroutine call to compute
    //     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
    //     this code.
    //
    INTEGER i = 0;
    for (i = 1; i <= k; i = i + 1) {
        dsigma[i - 1] = Rlamc3(dsigma[i - 1], dsigma[i - 1]) - dsigma[i - 1];
    }
    //
    //     Book keeping.
    //
    INTEGER iwk1 = 1;
    INTEGER iwk2 = iwk1 + k;
    INTEGER iwk3 = iwk2 + k;
    INTEGER iwk2i = iwk2 - 1;
    INTEGER iwk3i = iwk3 - 1;
    //
    //     Normalize Z.
    //
    REAL rho = Rnrm2(k, z, 1);
    Rlascl("G", 0, 0, rho, one, k, 1, z, k, info);
    rho = rho * rho;
    //
    //     Initialize WORK(IWK3).
    //
    Rlaset("A", k, 1, one, one, &work[iwk3 - 1], k);
    //
    //     Compute the updated singular values, the arrays DIFL, DIFR,
    //     and the updated Z.
    //
    INTEGER j = 0;
    for (j = 1; j <= k; j = j + 1) {
        Rlasd4(k, j, dsigma, z, &work[iwk1 - 1], rho, d[j - 1], &work[iwk2 - 1], info);
        //
        //        If the root finder fails, report the convergence failure.
        //
        if (info != 0) {
            return;
        }
        work[(iwk3i + j) - 1] = work[(iwk3i + j) - 1] * work[j - 1] * work[(iwk2i + j) - 1];
        difl[j - 1] = -work[j - 1];
        difr[(j - 1)] = -work[(j + 1) - 1];
        for (i = 1; i <= j - 1; i = i + 1) {
            work[(iwk3i + i) - 1] = work[(iwk3i + i) - 1] * work[i - 1] * work[(iwk2i + i) - 1] / (dsigma[i - 1] - dsigma[j - 1]) / (dsigma[i - 1] + dsigma[j - 1]);
        }
        for (i = j + 1; i <= k; i = i + 1) {
            work[(iwk3i + i) - 1] = work[(iwk3i + i) - 1] * work[i - 1] * work[(iwk2i + i) - 1] / (dsigma[i - 1] - dsigma[j - 1]) / (dsigma[i - 1] + dsigma[j - 1]);
        }
    }
    //
    //     Compute updated Z.
    //
    for (i = 1; i <= k; i = i + 1) {
        z[i - 1] = sign(sqrt(abs(work[(iwk3i + i) - 1])), z[i - 1]);
    }
    //
    //     Update VF and VL.
    //
    REAL diflj = 0.0;
    REAL dj = 0.0;
    REAL dsigj = 0.0;
    REAL difrj = 0.0;
    REAL dsigjp = 0.0;
    REAL temp = 0.0;
    for (j = 1; j <= k; j = j + 1) {
        diflj = difl[j - 1];
        dj = d[j - 1];
        dsigj = -dsigma[j - 1];
        if (j < k) {
            difrj = -difr[(j - 1)];
            dsigjp = -dsigma[(j + 1) - 1];
        }
        work[j - 1] = -z[j - 1] / diflj / (dsigma[j - 1] + dj);
        for (i = 1; i <= j - 1; i = i + 1) {
            work[i - 1] = z[i - 1] / Rlamc3(dsigma[i - 1], dsigma[i - 1]);
        }
        for (i = j + 1; i <= k; i = i + 1) {
            work[i - 1] = z[i - 1] / (Rlamc3(dsigma[i - 1], dsigjp) + difrj) / (dsigma[i - 1] + dj);
        }
        temp = Rnrm2(k, work, 1);
        work[(iwk2i + j) - 1] = Rdot(k, work, 1, vf, 1) / temp;
        work[(iwk3i + j) - 1] = Rdot(k, work, 1, vl, 1) / temp;
        if (icompq == 1) {
            difr[(j - 1) + (2 - 1) * lddifr] = temp;
        }
    }
    //
    Rcopy(k, &work[iwk2 - 1], 1, vf, 1);
    Rcopy(k, &work[iwk3 - 1], 1, vl, 1);
    //
    //     End of Rlasd8
    //
}

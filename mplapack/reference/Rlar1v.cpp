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

void Rlar1v(INTEGER const n, INTEGER const b1, INTEGER const bn, REAL const lambda, REAL *d, REAL *l, REAL *ld, REAL *lld, REAL const pivmin, REAL const gaptol, REAL *z, bool const wantnc, INTEGER &negcnt, REAL &ztz, REAL &mingma, INTEGER &r, INTEGER *isuppz, REAL &nrminv, REAL &resid, REAL &rqcorr, REAL *work) {
    REAL eps = 0.0;
    INTEGER r1 = 0;
    INTEGER r2 = 0;
    INTEGER indlpl = 0;
    INTEGER indumn = 0;
    INTEGER inds = 0;
    INTEGER indp = 0;
    const REAL zero = 0.0;
    bool sawnan1 = false;
    INTEGER neg1 = 0;
    REAL s = 0.0;
    INTEGER i = 0;
    REAL dplus = 0.0;
    bool sawnan2 = false;
    INTEGER neg2 = 0;
    REAL dminus = 0.0;
    REAL tmp = 0.0;
    const REAL one = 1.0;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    eps = Rlamch("Precision");
    //
    if (r == 0) {
        r1 = b1;
        r2 = bn;
    } else {
        r1 = r;
        r2 = r;
    }
    //
    //     Storage for LPLUS
    indlpl = 0;
    //     Storage for UMINUS
    indumn = n;
    inds = 2 * n + 1;
    indp = 3 * n + 1;
    //
    if (b1 == 1) {
        work[inds - 1] = zero;
    } else {
        work[(inds + b1 - 1) - 1] = lld[(b1 - 1) - 1];
    }
    //
    //     Compute the stationary transform (using the differential form)
    //     until the index R2.
    //
    sawnan1 = false;
    neg1 = 0;
    s = work[(inds + b1 - 1) - 1] - lambda;
    for (i = b1; i <= r1 - 1; i = i + 1) {
        dplus = d[i - 1] + s;
        work[(indlpl + i) - 1] = ld[i - 1] / dplus;
        if (dplus < zero) {
            neg1++;
        }
        work[(inds + i) - 1] = s * work[(indlpl + i) - 1] * l[i - 1];
        s = work[(inds + i) - 1] - lambda;
    }
    sawnan1 = Risnan(s);
    if (sawnan1) {
        goto statement_60;
    }
    for (i = r1; i <= r2 - 1; i = i + 1) {
        dplus = d[i - 1] + s;
        work[(indlpl + i) - 1] = ld[i - 1] / dplus;
        work[(inds + i) - 1] = s * work[(indlpl + i) - 1] * l[i - 1];
        s = work[(inds + i) - 1] - lambda;
    }
    sawnan1 = Risnan(s);
//
statement_60:
    if (sawnan1) {
        //        Runs a slower version of the above loop if a NaN is detected
        neg1 = 0;
        s = work[(inds + b1 - 1) - 1] - lambda;
        for (i = b1; i <= r1 - 1; i = i + 1) {
            dplus = d[i - 1] + s;
            if (abs(dplus) < pivmin) {
                dplus = -pivmin;
            }
            work[(indlpl + i) - 1] = ld[i - 1] / dplus;
            if (dplus < zero) {
                neg1++;
            }
            work[(inds + i) - 1] = s * work[(indlpl + i) - 1] * l[i - 1];
            if (work[(indlpl + i) - 1] == zero) {
                work[(inds + i) - 1] = lld[i - 1];
            }
            s = work[(inds + i) - 1] - lambda;
        }
        for (i = r1; i <= r2 - 1; i = i + 1) {
            dplus = d[i - 1] + s;
            if (abs(dplus) < pivmin) {
                dplus = -pivmin;
            }
            work[(indlpl + i) - 1] = ld[i - 1] / dplus;
            work[(inds + i) - 1] = s * work[(indlpl + i) - 1] * l[i - 1];
            if (work[(indlpl + i) - 1] == zero) {
                work[(inds + i) - 1] = lld[i - 1];
            }
            s = work[(inds + i) - 1] - lambda;
        }
    }
    //
    //     Compute the progressive transform (using the differential form)
    //     until the index R1
    //
    sawnan2 = false;
    neg2 = 0;
    work[(indp + bn - 1) - 1] = d[bn - 1] - lambda;
    for (i = bn - 1; i >= r1; i = i - 1) {
        dminus = lld[i - 1] + work[(indp + i) - 1];
        tmp = d[i - 1] / dminus;
        if (dminus < zero) {
            neg2++;
        }
        work[(indumn + i) - 1] = l[i - 1] * tmp;
        work[(indp + i - 1) - 1] = work[(indp + i) - 1] * tmp - lambda;
    }
    tmp = work[(indp + r1 - 1) - 1];
    sawnan2 = Risnan(tmp);
    //
    if (sawnan2) {
        //        Runs a slower version of the above loop if a NaN is detected
        neg2 = 0;
        for (i = bn - 1; i >= r1; i = i - 1) {
            dminus = lld[i - 1] + work[(indp + i) - 1];
            if (abs(dminus) < pivmin) {
                dminus = -pivmin;
            }
            tmp = d[i - 1] / dminus;
            if (dminus < zero) {
                neg2++;
            }
            work[(indumn + i) - 1] = l[i - 1] * tmp;
            work[(indp + i - 1) - 1] = work[(indp + i) - 1] * tmp - lambda;
            if (tmp == zero) {
                work[(indp + i - 1) - 1] = d[i - 1] - lambda;
            }
        }
    }
    //
    //     Find the index (from R1 to R2) of the largest (in magnitude)
    //     diagonal element of the inverse
    //
    mingma = work[(inds + r1 - 1) - 1] + work[(indp + r1 - 1) - 1];
    if (mingma < zero) {
        neg1++;
    }
    if (wantnc) {
        negcnt = neg1 + neg2;
    } else {
        negcnt = -1;
    }
    if (abs(mingma) == zero) {
        mingma = eps * work[(inds + r1 - 1) - 1];
    }
    r = r1;
    for (i = r1; i <= r2 - 1; i = i + 1) {
        tmp = work[(inds + i) - 1] + work[(indp + i) - 1];
        if (tmp == zero) {
            tmp = eps * work[(inds + i) - 1];
        }
        if (abs(tmp) <= abs(mingma)) {
            mingma = tmp;
            r = i + 1;
        }
    }
    //
    //     Compute the FP vector: solve N^T v = e_r
    //
    isuppz[1 - 1] = b1;
    isuppz[2 - 1] = bn;
    z[r - 1] = one;
    ztz = one;
    //
    //     Compute the FP vector upwards from R
    //
    if (!sawnan1 && !sawnan2) {
        for (i = r - 1; i >= b1; i = i - 1) {
            z[i - 1] = -(work[(indlpl + i) - 1] * z[(i + 1) - 1]);
            if ((abs(z[i - 1]) + abs(z[(i + 1) - 1])) * abs(ld[i - 1]) < gaptol) {
                z[i - 1] = zero;
                isuppz[1 - 1] = i + 1;
                goto statement_220;
            }
            ztz += z[i - 1] * z[i - 1];
        }
    statement_220:;
    } else {
        //        Run slower loop if NaN occurred.
        for (i = r - 1; i >= b1; i = i - 1) {
            if (z[(i + 1) - 1] == zero) {
                z[i - 1] = -(ld[(i + 1) - 1] / ld[i - 1]) * z[(i + 2) - 1];
            } else {
                z[i - 1] = -(work[(indlpl + i) - 1] * z[(i + 1) - 1]);
            }
            if ((abs(z[i - 1]) + abs(z[(i + 1) - 1])) * abs(ld[i - 1]) < gaptol) {
                z[i - 1] = zero;
                isuppz[1 - 1] = i + 1;
                goto statement_240;
            }
            ztz += z[i - 1] * z[i - 1];
        }
    statement_240:;
    }
    //
    //     Compute the FP vector downwards from R in blocks of size BLKSIZ
    if (!sawnan1 && !sawnan2) {
        for (i = r; i <= bn - 1; i = i + 1) {
            z[(i + 1) - 1] = -(work[(indumn + i) - 1] * z[i - 1]);
            if ((abs(z[i - 1]) + abs(z[(i + 1) - 1])) * abs(ld[i - 1]) < gaptol) {
                z[(i + 1) - 1] = zero;
                isuppz[2 - 1] = i;
                goto statement_260;
            }
            ztz += z[(i + 1) - 1] * z[(i + 1) - 1];
        }
    statement_260:;
    } else {
        //        Run slower loop if NaN occurred.
        for (i = r; i <= bn - 1; i = i + 1) {
            if (z[i - 1] == zero) {
                z[(i + 1) - 1] = -(ld[(i - 1) - 1] / ld[i - 1]) * z[(i - 1) - 1];
            } else {
                z[(i + 1) - 1] = -(work[(indumn + i) - 1] * z[i - 1]);
            }
            if ((abs(z[i - 1]) + abs(z[(i + 1) - 1])) * abs(ld[i - 1]) < gaptol) {
                z[(i + 1) - 1] = zero;
                isuppz[2 - 1] = i;
                goto statement_280;
            }
            ztz += z[(i + 1) - 1] * z[(i + 1) - 1];
        }
    statement_280:;
    }
    //
    //     Compute quantities for convergence test
    //
    tmp = one / ztz;
    nrminv = sqrt(tmp);
    resid = abs(mingma) * nrminv;
    rqcorr = mingma * tmp;
    //
    //     End of Rlar1v
    //
}

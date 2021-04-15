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

void Rggbal(const char *job, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, INTEGER &ilo, INTEGER &ihi, REAL *lscale, REAL *rscale, REAL *work, INTEGER &info) {
    const REAL one = 1.0;
    INTEGER i = 0;
    INTEGER k = 0;
    INTEGER l = 0;
    INTEGER lm1 = 0;
    INTEGER j = 0;
    INTEGER jp1 = 0;
    const REAL zero = 0.0;
    INTEGER m = 0;
    INTEGER iflow = 0;
    INTEGER ip1 = 0;
    INTEGER nr = 0;
    const REAL sclfac = 1.0e+1;
    REAL basl = 0.0;
    REAL tb = 0.0;
    REAL ta = 0.0;
    REAL coef = 0.0;
    REAL coef2 = 0.0;
    const REAL half = 0.5e+0;
    REAL coef5 = 0.0;
    INTEGER nrp2 = 0;
    REAL beta = 0.0;
    INTEGER it = 0;
    REAL gamma = 0.0;
    REAL ew = 0.0;
    REAL ewc = 0.0;
    REAL pgamma = 0.0;
    const REAL three = 3.0e+0;
    REAL t = 0.0;
    REAL tc = 0.0;
    INTEGER kount = 0;
    REAL sum = 0.0;
    REAL alpha = 0.0;
    REAL cmax = 0.0;
    REAL cor = 0.0;
    REAL sfmin = 0.0;
    REAL sfmax = 0.0;
    INTEGER lsfmin = 0;
    INTEGER lsfmax = 0;
    INTEGER irab = 0;
    REAL rab = 0.0;
    INTEGER lrab = 0;
    INTEGER ir = 0;
    INTEGER icab = 0;
    REAL cab = 0.0;
    INTEGER lcab = 0;
    INTEGER jc = 0;
    //
    //  -- LAPACK computational routine --
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
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    if (!Mlsame(job, "N") && !Mlsame(job, "P") && !Mlsame(job, "S") && !Mlsame(job, "B")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rggbal", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        ilo = 1;
        ihi = n;
        return;
    }
    //
    if (n == 1) {
        ilo = 1;
        ihi = n;
        lscale[1 - 1] = one;
        rscale[1 - 1] = one;
        return;
    }
    //
    if (Mlsame(job, "N")) {
        ilo = 1;
        ihi = n;
        for (i = 1; i <= n; i = i + 1) {
            lscale[i - 1] = one;
            rscale[i - 1] = one;
        }
        return;
    }
    //
    k = 1;
    l = n;
    if (Mlsame(job, "S")) {
        goto statement_190;
    }
    //
    goto statement_30;
//
//     Permute the matrices A and B to isolate the eigenvalues.
//
//     Find row with one nonzero in columns 1 through L
//
statement_20:
    l = lm1;
    if (l != 1) {
        goto statement_30;
    }
    //
    rscale[1 - 1] = one;
    lscale[1 - 1] = one;
    goto statement_190;
//
statement_30:
    lm1 = l - 1;
    for (i = l; i >= 1; i = i - 1) {
        for (j = 1; j <= lm1; j = j + 1) {
            jp1 = j + 1;
            if (a[(i - 1) + (j - 1) * lda] != zero || b[(i - 1) + (j - 1) * ldb] != zero) {
                goto statement_50;
            }
        }
        j = l;
        goto statement_70;
    //
    statement_50:
        for (j = jp1; j <= l; j = j + 1) {
            if (a[(i - 1) + (j - 1) * lda] != zero || b[(i - 1) + (j - 1) * ldb] != zero) {
                goto statement_80;
            }
        }
        j = jp1 - 1;
    //
    statement_70:
        m = l;
        iflow = 1;
        goto statement_160;
    statement_80:;
    }
    goto statement_100;
//
//     Find column with one nonzero in rows K through N
//
statement_90:
    k++;
//
statement_100:
    for (j = k; j <= l; j = j + 1) {
        for (i = k; i <= lm1; i = i + 1) {
            ip1 = i + 1;
            if (a[(i - 1) + (j - 1) * lda] != zero || b[(i - 1) + (j - 1) * ldb] != zero) {
                goto statement_120;
            }
        }
        i = l;
        goto statement_140;
    statement_120:
        for (i = ip1; i <= l; i = i + 1) {
            if (a[(i - 1) + (j - 1) * lda] != zero || b[(i - 1) + (j - 1) * ldb] != zero) {
                goto statement_150;
            }
        }
        i = ip1 - 1;
    statement_140:
        m = k;
        iflow = 2;
        goto statement_160;
    statement_150:;
    }
    goto statement_190;
//
//     Permute rows M and I
//
statement_160:
    lscale[m - 1] = i;
    if (i == m) {
        goto statement_170;
    }
    Rswap(n - k + 1, &a[(i - 1) + (k - 1) * lda], lda, &a[(m - 1) + (k - 1) * lda], lda);
    Rswap(n - k + 1, &b[(i - 1) + (k - 1) * ldb], ldb, &b[(m - 1) + (k - 1) * ldb], ldb);
//
//     Permute columns M and J
//
statement_170:
    rscale[m - 1] = j;
    if (j == m) {
        goto statement_180;
    }
    Rswap(l, &a[(j - 1) * lda], 1, &a[(m - 1) * lda], 1);
    Rswap(l, &b[(j - 1) * ldb], 1, &b[(m - 1) * ldb], 1);
//
statement_180:
    switch (iflow) {
    case 1:
        goto statement_20;
    case 2:
        goto statement_90;
    default:
        break;
    }
//
statement_190:
    ilo = k;
    ihi = l;
    //
    if (Mlsame(job, "P")) {
        for (i = ilo; i <= ihi; i = i + 1) {
            lscale[i - 1] = one;
            rscale[i - 1] = one;
        }
        return;
    }
    //
    if (ilo == ihi) {
        return;
    }
    //
    //     Balance the submatrix in rows ILO to IHI.
    //
    nr = ihi - ilo + 1;
    for (i = ilo; i <= ihi; i = i + 1) {
        rscale[i - 1] = zero;
        lscale[i - 1] = zero;
        //
        work[i - 1] = zero;
        work[(i + n) - 1] = zero;
        work[(i + 2 * n) - 1] = zero;
        work[(i + 3 * n) - 1] = zero;
        work[(i + 4 * n) - 1] = zero;
        work[(i + 5 * n) - 1] = zero;
    }
    //
    //     Compute right side vector in resulting linear equations
    //
    basl = log10(sclfac);
    for (i = ilo; i <= ihi; i = i + 1) {
        for (j = ilo; j <= ihi; j = j + 1) {
            tb = b[(i - 1) + (j - 1) * ldb];
            ta = a[(i - 1) + (j - 1) * lda];
            if (ta == zero) {
                goto statement_210;
            }
            ta = log10(abs(ta)) / basl;
        statement_210:
            if (tb == zero) {
                goto statement_220;
            }
            tb = log10(abs(tb)) / basl;
        statement_220:
            work[(i + 4 * n) - 1] = work[(i + 4 * n) - 1] - ta - tb;
            work[(j + 5 * n) - 1] = work[(j + 5 * n) - 1] - ta - tb;
        }
    }
    //
    coef = one / castREAL(2 * nr);
    coef2 = coef * coef;
    coef5 = half * coef2;
    nrp2 = nr + 2;
    beta = zero;
    it = 1;
//
//     Start generalized conjugate gradient iteration
//
statement_250:
    //
    gamma = Rdot(nr, &work[(ilo + 4 * n) - 1], 1, &work[(ilo + 4 * n) - 1], 1) + Rdot(nr, &work[(ilo + 5 * n) - 1], 1, &work[(ilo + 5 * n) - 1], 1);
    //
    ew = zero;
    ewc = zero;
    for (i = ilo; i <= ihi; i = i + 1) {
        ew += work[(i + 4 * n) - 1];
        ewc += work[(i + 5 * n) - 1];
    }
    //
    gamma = coef * gamma - coef2 * (pow2(ew) + pow2(ewc)) - coef5 * pow2((ew - ewc));
    if (gamma == zero) {
        goto statement_350;
    }
    if (it != 1) {
        beta = gamma / pgamma;
    }
    t = coef5 * (ewc - three * ew);
    tc = coef5 * (ew - three * ewc);
    //
    Rscal(nr, beta, &work[ilo - 1], 1);
    Rscal(nr, beta, &work[(ilo + n) - 1], 1);
    //
    Raxpy(nr, coef, &work[(ilo + 4 * n) - 1], 1, &work[(ilo + n) - 1], 1);
    Raxpy(nr, coef, &work[(ilo + 5 * n) - 1], 1, &work[ilo - 1], 1);
    //
    for (i = ilo; i <= ihi; i = i + 1) {
        work[i - 1] += tc;
        work[(i + n) - 1] += t;
    }
    //
    //     Apply matrix to vector
    //
    for (i = ilo; i <= ihi; i = i + 1) {
        kount = 0;
        sum = zero;
        for (j = ilo; j <= ihi; j = j + 1) {
            if (a[(i - 1) + (j - 1) * lda] == zero) {
                goto statement_280;
            }
            kount++;
            sum += work[j - 1];
        statement_280:
            if (b[(i - 1) + (j - 1) * ldb] == zero) {
                goto statement_290;
            }
            kount++;
            sum += work[j - 1];
        statement_290:;
        }
        work[(i + 2 * n) - 1] = castREAL(kount) * work[(i + n) - 1] + sum;
    }
    //
    for (j = ilo; j <= ihi; j = j + 1) {
        kount = 0;
        sum = zero;
        for (i = ilo; i <= ihi; i = i + 1) {
            if (a[(i - 1) + (j - 1) * lda] == zero) {
                goto statement_310;
            }
            kount++;
            sum += work[(i + n) - 1];
        statement_310:
            if (b[(i - 1) + (j - 1) * ldb] == zero) {
                goto statement_320;
            }
            kount++;
            sum += work[(i + n) - 1];
        statement_320:;
        }
        work[(j + 3 * n) - 1] = castREAL(kount) * work[j - 1] + sum;
    }
    //
    sum = Rdot(nr, &work[(ilo + n) - 1], 1, &work[(ilo + 2 * n) - 1], 1) + Rdot(nr, &work[ilo - 1], 1, &work[(ilo + 3 * n) - 1], 1);
    alpha = gamma / sum;
    //
    //     Determine correction to current iteration
    //
    cmax = zero;
    for (i = ilo; i <= ihi; i = i + 1) {
        cor = alpha * work[(i + n) - 1];
        if (abs(cor) > cmax) {
            cmax = abs(cor);
        }
        lscale[i - 1] += cor;
        cor = alpha * work[i - 1];
        if (abs(cor) > cmax) {
            cmax = abs(cor);
        }
        rscale[i - 1] += cor;
    }
    if (cmax < half) {
        goto statement_350;
    }
    //
    Raxpy(nr, -alpha, &work[(ilo + 2 * n) - 1], 1, &work[(ilo + 4 * n) - 1], 1);
    Raxpy(nr, -alpha, &work[(ilo + 3 * n) - 1], 1, &work[(ilo + 5 * n) - 1], 1);
    //
    pgamma = gamma;
    it++;
    if (it <= nrp2) {
        goto statement_250;
    }
//
//     End generalized conjugate gradient iteration
//
statement_350:
    sfmin = Rlamch("S");
    sfmax = one / sfmin;
    lsfmin = castINTEGER(log10(sfmin - 1) / basl + one);
    lsfmax = castINTEGER(log10(sfmax - 1) / basl);
    for (i = ilo; i <= ihi; i = i + 1) {
        irab = iRamax(n - ilo + 1, &a[(i - 1) + (ilo - 1) * lda], lda);
        rab = abs(a[(i - 1) + ((irab + ilo - 1) - 1) * lda]);
        irab = iRamax(n - ilo + 1, &b[(i - 1) + (ilo - 1) * ldb], ldb);
        rab = max(rab, abs(b[(i - 1) + ((irab + ilo - 1) - 1) * ldb]));
        lrab = castINTEGER(log10(rab + sfmin) / basl + one);
        ir = castINTEGER(lscale[i - 1] + sign(half, lscale[i - 1]));
        ir = min({max(ir, lsfmin), lsfmax, lsfmax - lrab});
        lscale[i - 1] = pow(sclfac, ir);
        icab = iRamax(ihi, &a[(i - 1) * lda], 1);
        cab = abs(a[(icab - 1) + (i - 1) * lda]);
        icab = iRamax(ihi, &b[(i - 1) * ldb], 1);
        cab = max(cab, abs(b[(icab - 1) + (i - 1) * ldb]));
        lcab = castINTEGER(log10((cab + sfmin)) / basl + one);
        jc = castINTEGER(rscale[i - 1] + sign(half, rscale[i - 1]));
        jc = min({max(jc, lsfmin), lsfmax, lsfmax - lcab});
        rscale[i - 1] = pow(sclfac, jc);
    }
    //
    //     Row scaling of matrices A and B
    //
    for (i = ilo; i <= ihi; i = i + 1) {
        Rscal(n - ilo + 1, lscale[i - 1], &a[(i - 1) + (ilo - 1) * lda], lda);
        Rscal(n - ilo + 1, lscale[i - 1], &b[(i - 1) + (ilo - 1) * ldb], ldb);
    }
    //
    //     Column scaling of matrices A and B
    //
    for (j = ilo; j <= ihi; j = j + 1) {
        Rscal(ihi, rscale[j - 1], &a[(j - 1) * lda], 1);
        Rscal(ihi, rscale[j - 1], &b[(j - 1) * ldb], 1);
    }
    //
    //     End of Rggbal
    //
}

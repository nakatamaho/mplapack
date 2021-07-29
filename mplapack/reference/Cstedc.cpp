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

void Cstedc(const char *compz, INTEGER const n, REAL *d, REAL *e, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
    bool lquery = false;
    INTEGER icompz = 0;
    INTEGER smlsiz = 0;
    INTEGER lwmin = 0;
    INTEGER liwmin = 0;
    INTEGER lrwmin = 0;
    const REAL two = 2.0;
    INTEGER lgn = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    INTEGER ll = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    REAL orgnrm = 0.0;
    REAL eps = 0.0;
    INTEGER start = 0;
    INTEGER finish = 0;
    REAL tiny = 0.0;
    INTEGER m = 0;
    INTEGER ii = 0;
    INTEGER k = 0;
    REAL p = 0.0;
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
    //     Test the input parameters.
    //
    info = 0;
    lquery = (lwork == -1 || lrwork == -1 || liwork == -1);
    //
    if (Mlsame(compz, "N")) {
        icompz = 0;
    } else if (Mlsame(compz, "V")) {
        icompz = 1;
    } else if (Mlsame(compz, "I")) {
        icompz = 2;
    } else {
        icompz = -1;
    }
    if (icompz < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if ((ldz < 1) || (icompz > 0 && ldz < max((INTEGER)1, n))) {
        info = -6;
    }
    //
    if (info == 0) {
        //
        //        Compute the workspace requirements
        //
        smlsiz = iMlaenv(9, "Cstedc", " ", 0, 0, 0, 0);
        if (n <= 1 || icompz == 0) {
            lwmin = 1;
            liwmin = 1;
            lrwmin = 1;
        } else if (n <= smlsiz) {
            lwmin = 1;
            liwmin = 1;
            lrwmin = 2 * (n - 1);
        } else if (icompz == 1) {
            lgn = castINTEGER(log(castREAL(n)) / log(two));
            if (pow(2, lgn) < n) {
                lgn++;
            }
            if (pow(2, lgn) < n) {
                lgn++;
            }
            lwmin = n * n;
            lrwmin = 1 + 3 * n + 2 * n * lgn + 4 * n * n;
            liwmin = 6 + 6 * n + 5 * n * lgn;
        } else if (icompz == 2) {
            lwmin = 1;
            lrwmin = 1 + 4 * n + 2 * n * n;
            liwmin = 3 + 5 * n;
        }
        work[1 - 1] = lwmin;
        rwork[1 - 1] = lrwmin;
        iwork[1 - 1] = liwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -8;
        } else if (lrwork < lrwmin && !lquery) {
            info = -10;
        } else if (liwork < liwmin && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cstedc", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    if (n == 1) {
        if (icompz != 0) {
            z[(1 - 1)] = one;
        }
        return;
    }
    //
    //     If the following conditional clause is removed, then the routine
    //     will use the Divide and Conquer routine to compute only the
    //     eigenvalues, which requires (3N + 3N**2) real workspace and
    //     (2 + 5N + 2N lg(N)) integer workspace.
    //     Since on many architectures Rsterf is much faster than any other
    //     algorithm for finding eigenvalues only, it is used here
    //     as the default. If the conditional clause is removed, then
    //     information on the size of workspace needs to be changed.
    //
    //     If COMPZ = 'N', use Rsterf to compute the eigenvalues.
    //
    if (icompz == 0) {
        Rsterf(n, d, e, info);
        goto statement_70;
    }
    //
    //     If N is smaller than the minimum divide size (SMLSIZ+1), then
    //     solve the problem with another solver.
    //
    if (n <= smlsiz) {
        //
        Csteqr(compz, n, d, e, z, ldz, rwork, info);
        //
    } else {
        //
        //        If COMPZ = 'I', we simply call Rstedc instead.
        //
        if (icompz == 2) {
            Rlaset("Full", n, n, zero, one, rwork, n);
            ll = n * n + 1;
            Rstedc("I", n, d, e, rwork, n, &rwork[ll - 1], lrwork - ll + 1, iwork, liwork, info);
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    z[(i - 1) + (j - 1) * ldz] = rwork[((j - 1) * n + i) - 1];
                }
            }
            goto statement_70;
        }
        //
        //        From now on, only option left to be handled is COMPZ = 'V',
        //        i.e. ICOMPZ = 1.
        //
        //        Scale.
        //
        orgnrm = Rlanst("M", n, d, e);
        if (orgnrm == zero) {
            goto statement_70;
        }
        //
        eps = Rlamch("Epsilon");
        //
        start = 1;
    //
    //        while ( START <= N )
    //
    statement_30:
        if (start <= n) {
            //
            //           Let FINISH be the position of the next subdiagonal entry
            //           such that E( FINISH ) <= TINY or FINISH = N if no such
            //           subdiagonal exists.  The matrix identified by the elements
            //           between START and FINISH constitutes an independent
            //           sub-problem.
            //
            finish = start;
        statement_40:
            if (finish < n) {
                tiny = eps * sqrt(abs(d[finish - 1])) * sqrt(abs(d[(finish + 1) - 1]));
                if (abs(e[finish - 1]) > tiny) {
                    finish++;
                    goto statement_40;
                }
            }
            //
            //           (Sub) Problem determined.  Compute its size and solve it.
            //
            m = finish - start + 1;
            if (m > smlsiz) {
                //
                //              Scale.
                //
                orgnrm = Rlanst("M", m, &d[start - 1], &e[start - 1]);
                Rlascl("G", 0, 0, orgnrm, one, m, 1, &d[start - 1], m, info);
                Rlascl("G", 0, 0, orgnrm, one, m - 1, 1, &e[start - 1], m - 1, info);
                //
                Claed0(n, m, &d[start - 1], &e[start - 1], &z[(start - 1) * ldz], ldz, work, n, rwork, iwork, info);
                if (info > 0) {
                    info = (info / (m + 1) + start - 1) * (n + 1) + mod(info, (m + 1)) + start - 1;
                    goto statement_70;
                }
                //
                //              Scale back.
                //
                Rlascl("G", 0, 0, one, orgnrm, m, 1, &d[start - 1], m, info);
                //
            } else {
                Rsteqr("I", m, &d[start - 1], &e[start - 1], rwork, m, &rwork[(m * m + 1) - 1], info);
                Clacrm(n, m, &z[(start - 1) * ldz], ldz, rwork, m, work, n, &rwork[(m * m + 1) - 1]);
                Clacpy("A", n, m, work, n, &z[(start - 1) * ldz], ldz);
                if (info > 0) {
                    info = start * (n + 1) + finish;
                    goto statement_70;
                }
            }
            //
            start = finish + 1;
            goto statement_30;
        }
        //
        //        endwhile
        //
        //        Use Selection Sort to minimize swaps of eigenvectors
        //
        for (ii = 2; ii <= n; ii = ii + 1) {
            i = ii - 1;
            k = i;
            p = d[i - 1];
            for (j = ii; j <= n; j = j + 1) {
                if (d[j - 1] < p) {
                    k = j;
                    p = d[j - 1];
                }
            }
            if (k != i) {
                d[k - 1] = d[i - 1];
                d[i - 1] = p;
                Cswap(n, &z[(i - 1) * ldz], 1, &z[(k - 1) * ldz], 1);
            }
        }
    }
//
statement_70:
    work[1 - 1] = lwmin;
    rwork[1 - 1] = lrwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Cstedc
    //
}

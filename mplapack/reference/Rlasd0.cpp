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

void Rlasd0(INTEGER const n, INTEGER const sqre, REAL *d, REAL *e, REAL *u, INTEGER const ldu, REAL *vt, INTEGER const ldvt, INTEGER const smlsiz, INTEGER *iwork, REAL *work, INTEGER &info) {
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
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    //
    if (n < 0) {
        info = -1;
    } else if ((sqre < 0) || (sqre > 1)) {
        info = -2;
    }
    //
    INTEGER m = n + sqre;
    //
    if (ldu < n) {
        info = -6;
    } else if (ldvt < m) {
        info = -8;
    } else if (smlsiz < 3) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Rlasd0", -info);
        return;
    }
    //
    //     If the input matrix is too small, call Rlasdq to find the SVD.
    //
    if (n <= smlsiz) {
        Rlasdq("U", sqre, n, m, n, 0, d, e, vt, ldvt, u, ldu, u, ldu, work, info);
        return;
    }
    //
    //     Set up the computation tree.
    //
    INTEGER inode = 1;
    INTEGER ndiml = inode + n;
    INTEGER ndimr = ndiml + n;
    INTEGER idxq = ndimr + n;
    INTEGER iwk = idxq + n;
    INTEGER nlvl = 0;
    INTEGER nd = 0;
    Rlasdt(n, nlvl, nd, &iwork[inode - 1], &iwork[ndiml - 1], &iwork[ndimr - 1], smlsiz);
    //
    //     For the nodes on bottom level of the tree, solve
    //     their subproblems by Rlasdq.
    //
    INTEGER ndb1 = (nd + 1) / 2;
    INTEGER ncc = 0;
    INTEGER i = 0;
    INTEGER i1 = 0;
    INTEGER ic = 0;
    INTEGER nl = 0;
    INTEGER nlp1 = 0;
    INTEGER nr = 0;
    INTEGER nrp1 = 0;
    INTEGER nlf = 0;
    INTEGER nrf = 0;
    INTEGER sqrei = 0;
    INTEGER itemp = 0;
    INTEGER j = 0;
    for (i = ndb1; i <= nd; i = i + 1) {
        //
        //     IC : center row of each node
        //     NL : number of rows of left  subproblem
        //     NR : number of rows of right subproblem
        //     NLF: starting row of the left   subproblem
        //     NRF: starting row of the right  subproblem
        //
        i1 = i - 1;
        ic = iwork[(inode + i1) - 1];
        nl = iwork[(ndiml + i1) - 1];
        nlp1 = nl + 1;
        nr = iwork[(ndimr + i1) - 1];
        nrp1 = nr + 1;
        nlf = ic - nl;
        nrf = ic + 1;
        sqrei = 1;
        Rlasdq("U", sqrei, nl, nlp1, nl, ncc, &d[nlf - 1], &e[nlf - 1], &vt[(nlf - 1) + (nlf - 1) * ldvt], ldvt, &u[(nlf - 1) + (nlf - 1) * ldu], ldu, &u[(nlf - 1) + (nlf - 1) * ldu], ldu, work, info);
        if (info != 0) {
            return;
        }
        itemp = idxq + nlf - 2;
        for (j = 1; j <= nl; j = j + 1) {
            iwork[(itemp + j) - 1] = j;
        }
        if (i == nd) {
            sqrei = sqre;
        } else {
            sqrei = 1;
        }
        nrp1 = nr + sqrei;
        Rlasdq("U", sqrei, nr, nrp1, nr, ncc, &d[nrf - 1], &e[nrf - 1], &vt[(nrf - 1) + (nrf - 1) * ldvt], ldvt, &u[(nrf - 1) + (nrf - 1) * ldu], ldu, &u[(nrf - 1) + (nrf - 1) * ldu], ldu, work, info);
        if (info != 0) {
            return;
        }
        itemp = idxq + ic;
        for (j = 1; j <= nr; j = j + 1) {
            iwork[(itemp + j - 1) - 1] = j;
        }
    }
    //
    //     Now conquer each subproblem bottom-up.
    //
    INTEGER lvl = 0;
    INTEGER lf = 0;
    INTEGER ll = 0;
    INTEGER im1 = 0;
    INTEGER idxqc = 0;
    REAL alpha = 0.0;
    REAL beta = 0.0;
    for (lvl = nlvl; lvl >= 1; lvl = lvl - 1) {
        //
        //        Find the first node LF and last node LL on the
        //        current level LVL.
        //
        if (lvl == 1) {
            lf = 1;
            ll = 1;
        } else {
            lf = pow(2, (lvl - 1));
            ll = 2 * lf - 1;
        }
        for (i = lf; i <= ll; i = i + 1) {
            im1 = i - 1;
            ic = iwork[(inode + im1) - 1];
            nl = iwork[(ndiml + im1) - 1];
            nr = iwork[(ndimr + im1) - 1];
            nlf = ic - nl;
            if ((sqre == 0) && (i == ll)) {
                sqrei = sqre;
            } else {
                sqrei = 1;
            }
            idxqc = idxq + nlf - 1;
            alpha = d[ic - 1];
            beta = e[ic - 1];
            Rlasd1(nl, nr, sqrei, &d[nlf - 1], alpha, beta, &u[(nlf - 1) + (nlf - 1) * ldu], ldu, &vt[(nlf - 1) + (nlf - 1) * ldvt], ldvt, &iwork[idxqc - 1], &iwork[iwk - 1], work, info);
            //
            //        Report the possible convergence failure.
            //
            if (info != 0) {
                return;
            }
        }
    }
    //
    //     End of Rlasd0
    //
}

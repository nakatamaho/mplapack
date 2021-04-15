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

void Rlasda(INTEGER const icompq, INTEGER const smlsiz, INTEGER const n, INTEGER const sqre, REAL *d, REAL *e, REAL *u, INTEGER const ldu, REAL *vt, INTEGER *k, REAL *difl, REAL *difr, REAL *z, REAL *poles, INTEGER *givptr, INTEGER *givcol, INTEGER const ldgcol, INTEGER *perm, REAL *givnum, REAL *c, REAL *s, REAL *work, INTEGER *iwork, INTEGER &info) {
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
    //     ..
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
    if ((icompq < 0) || (icompq > 1)) {
        info = -1;
    } else if (smlsiz < 3) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if ((sqre < 0) || (sqre > 1)) {
        info = -4;
    } else if (ldu < (n + sqre)) {
        info = -8;
    } else if (ldgcol < n) {
        info = -17;
    }
    if (info != 0) {
        Mxerbla("Rlasda", -info);
        return;
    }
    //
    INTEGER m = n + sqre;
    //
    //     If the input matrix is too small, call Rlasdq to find the SVD.
    //
    if (n <= smlsiz) {
        if (icompq == 0) {
            Rlasdq("U", sqre, n, 0, 0, 0, d, e, vt, ldu, u, ldu, u, ldu, work, info);
        } else {
            Rlasdq("U", sqre, n, m, n, 0, d, e, vt, ldu, u, ldu, u, ldu, work, info);
        }
        return;
    }
    //
    //     Book-keeping and  set up the computation tree.
    //
    INTEGER inode = 1;
    INTEGER ndiml = inode + n;
    INTEGER ndimr = ndiml + n;
    INTEGER idxq = ndimr + n;
    INTEGER iwk = idxq + n;
    //
    INTEGER ncc = 0;
    INTEGER nru = 0;
    //
    INTEGER smlszp = smlsiz + 1;
    INTEGER vf = 1;
    INTEGER vl = vf + m;
    INTEGER nwork1 = vl + m;
    INTEGER nwork2 = nwork1 + smlszp * smlszp;
    //
    INTEGER nlvl = 0;
    INTEGER nd = 0;
    Rlasdt(n, nlvl, nd, &iwork[inode - 1], &iwork[ndiml - 1], &iwork[ndimr - 1], smlsiz);
    //
    //     for the nodes on bottom level of the tree, solve
    //     their subproblems by Rlasdq.
    //
    INTEGER ndb1 = (nd + 1) / 2;
    INTEGER i = 0;
    INTEGER i1 = 0;
    INTEGER ic = 0;
    INTEGER nl = 0;
    INTEGER nlp1 = 0;
    INTEGER nr = 0;
    INTEGER nlf = 0;
    INTEGER nrf = 0;
    INTEGER idxqi = 0;
    INTEGER vfi = 0;
    INTEGER vli = 0;
    INTEGER sqrei = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER itemp = 0;
    INTEGER j = 0;
    INTEGER nrp1 = 0;
    for (i = ndb1; i <= nd; i = i + 1) {
        //
        //        IC : center row of each node
        //        NL : number of rows of left  subproblem
        //        NR : number of rows of right subproblem
        //        NLF: starting row of the left   subproblem
        //        NRF: starting row of the right  subproblem
        //
        i1 = i - 1;
        ic = iwork[(inode + i1) - 1];
        nl = iwork[(ndiml + i1) - 1];
        nlp1 = nl + 1;
        nr = iwork[(ndimr + i1) - 1];
        nlf = ic - nl;
        nrf = ic + 1;
        idxqi = idxq + nlf - 2;
        vfi = vf + nlf - 1;
        vli = vl + nlf - 1;
        sqrei = 1;
        if (icompq == 0) {
            Rlaset("A", nlp1, nlp1, zero, one, &work[nwork1 - 1], smlszp);
            Rlasdq("U", sqrei, nl, nlp1, nru, ncc, &d[nlf - 1], &e[nlf - 1], &work[nwork1 - 1], smlszp, &work[nwork2 - 1], nl, &work[nwork2 - 1], nl, &work[nwork2 - 1], info);
            itemp = nwork1 + nl * smlszp;
            Rcopy(nlp1, &work[nwork1 - 1], 1, &work[vfi - 1], 1);
            Rcopy(nlp1, &work[itemp - 1], 1, &work[vli - 1], 1);
        } else {
            Rlaset("A", nl, nl, zero, one, &u[(nlf - 1)], ldu);
            Rlaset("A", nlp1, nlp1, zero, one, &vt[(nlf - 1)], ldu);
            Rlasdq("U", sqrei, nl, nlp1, nl, ncc, &d[nlf - 1], &e[nlf - 1], &vt[(nlf - 1)], ldu, &u[(nlf - 1)], ldu, &u[(nlf - 1)], ldu, &work[nwork1 - 1], info);
            Rcopy(nlp1, &vt[(nlf - 1)], 1, &work[vfi - 1], 1);
            Rcopy(nlp1, &vt[(nlf - 1) + (nlp1 - 1) * ldvt], 1, &work[vli - 1], 1);
        }
        if (info != 0) {
            return;
        }
        for (j = 1; j <= nl; j = j + 1) {
            iwork[(idxqi + j) - 1] = j;
        }
        if ((i == nd) && (sqre == 0)) {
            sqrei = 0;
        } else {
            sqrei = 1;
        }
        idxqi += nlp1;
        vfi += nlp1;
        vli += nlp1;
        nrp1 = nr + sqrei;
        if (icompq == 0) {
            Rlaset("A", nrp1, nrp1, zero, one, &work[nwork1 - 1], smlszp);
            Rlasdq("U", sqrei, nr, nrp1, nru, ncc, &d[nrf - 1], &e[nrf - 1], &work[nwork1 - 1], smlszp, &work[nwork2 - 1], nr, &work[nwork2 - 1], nr, &work[nwork2 - 1], info);
            itemp = nwork1 + (nrp1 - 1) * smlszp;
            Rcopy(nrp1, &work[nwork1 - 1], 1, &work[vfi - 1], 1);
            Rcopy(nrp1, &work[itemp - 1], 1, &work[vli - 1], 1);
        } else {
            Rlaset("A", nr, nr, zero, one, &u[(nrf - 1)], ldu);
            Rlaset("A", nrp1, nrp1, zero, one, &vt[(nrf - 1)], ldu);
            Rlasdq("U", sqrei, nr, nrp1, nr, ncc, &d[nrf - 1], &e[nrf - 1], &vt[(nrf - 1)], ldu, &u[(nrf - 1)], ldu, &u[(nrf - 1)], ldu, &work[nwork1 - 1], info);
            Rcopy(nrp1, &vt[(nrf - 1)], 1, &work[vfi - 1], 1);
            Rcopy(nrp1, &vt[(nrf - 1) + (nrp1 - 1) * ldvt], 1, &work[vli - 1], 1);
        }
        if (info != 0) {
            return;
        }
        for (j = 1; j <= nr; j = j + 1) {
            iwork[(idxqi + j) - 1] = j;
        }
    }
    //
    //     Now conquer each subproblem bottom-up.
    //
    j = pow(2, nlvl);
    INTEGER lvl = 0;
    INTEGER lvl2 = 0;
    INTEGER lf = 0;
    INTEGER ll = 0;
    INTEGER im1 = 0;
    REAL alpha = 0.0;
    REAL beta = 0.0;
    for (lvl = nlvl; lvl >= 1; lvl = lvl - 1) {
        lvl2 = lvl * 2 - 1;
        //
        //        Find the first node LF and last node LL on
        //        the current level LVL.
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
            nrf = ic + 1;
            if (i == ll) {
                sqrei = sqre;
            } else {
                sqrei = 1;
            }
            vfi = vf + nlf - 1;
            vli = vl + nlf - 1;
            idxqi = idxq + nlf - 1;
            alpha = d[ic - 1];
            beta = e[ic - 1];
            if (icompq == 0) {
                Rlasd6(icompq, nl, nr, sqrei, &d[nlf - 1], &work[vfi - 1], &work[vli - 1], alpha, beta, iwork[idxqi - 1], perm, givptr[1 - 1], givcol, ldgcol, givnum, ldu, poles, difl, difr, z, k[1 - 1], &c[1 - 1], s[1 - 1], &work[nwork1 - 1], iwork[iwk - 1], info);
            } else {
                j = j - 1;
                Rlasd6(icompq, nl, nr, sqrei, &d[nlf - 1], &work[vfi - 1], &work[vli - 1], alpha, beta, iwork[idxqi - 1], perm[(nlf - 1) + (lvl - 1) * ldperm], givptr[j - 1], givcol[(nlf - 1) + (lvl2 - 1) * ldgivcol], ldgcol, givnum[(nlf - 1) + (lvl2 - 1) * ldgivnum], ldu, poles[(nlf - 1) + (lvl2 - 1) * ldpoles], difl[(nlf - 1) + (lvl - 1) * lddifl], difr[(nlf - 1) + (lvl2 - 1) * lddifr], &z[(nlf - 1) + (lvl - 1) * ldz], k[j - 1], &c[j - 1], s[j - 1], &work[nwork1 - 1], iwork[iwk - 1], info);
            }
            if (info != 0) {
                return;
            }
        }
    }
    //
    //     End of Rlasda
    //
}

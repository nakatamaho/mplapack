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

void Rlalsa(INTEGER const icompq, INTEGER const smlsiz, INTEGER const n, INTEGER const nrhs, REAL *b, INTEGER const ldb, REAL *bx, INTEGER const ldbx, REAL *u, INTEGER const ldu, REAL *vt, INTEGER *k, REAL *difl, REAL *difr, REAL *z, REAL *poles, INTEGER *givptr, INTEGER *givcol, INTEGER const ldgcol, INTEGER *perm, REAL *givnum, REAL *c, REAL *s, REAL *work, INTEGER *iwork, INTEGER &info) {
    INTEGER inode = 0;
    INTEGER ndiml = 0;
    INTEGER ndimr = 0;
    INTEGER nlvl = 0;
    INTEGER nd = 0;
    INTEGER ndb1 = 0;
    INTEGER i = 0;
    INTEGER i1 = 0;
    INTEGER ic = 0;
    INTEGER nl = 0;
    INTEGER nr = 0;
    INTEGER nlf = 0;
    INTEGER nrf = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER sqre = 0;
    INTEGER lvl = 0;
    INTEGER lvl2 = 0;
    INTEGER lf = 0;
    INTEGER ll = 0;
    INTEGER im1 = 0;
    INTEGER nlp1 = 0;
    INTEGER nrp1 = 0;
    INTEGER ldvt = ldu;
    INTEGER lddifl = ldu;
    INTEGER lddifr = ldu;
    INTEGER ldz = ldu;
    INTEGER ldpoles = ldu;
    INTEGER ldgivcol = ldgcol;
    INTEGER ldperm = ldgcol;
    INTEGER ldgivnum = ldu;
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
    } else if (n < smlsiz) {
        info = -3;
    } else if (nrhs < 1) {
        info = -4;
    } else if (ldb < n) {
        info = -6;
    } else if (ldbx < n) {
        info = -8;
    } else if (ldu < n) {
        info = -10;
    } else if (ldgcol < n) {
        info = -19;
    }
    if (info != 0) {
        Mxerbla("Rlalsa", -info);
        return;
    }
    //
    //     Book-keeping and  setting up the computation tree.
    //
    inode = 1;
    ndiml = inode + n;
    ndimr = ndiml + n;
    //
    Rlasdt(n, nlvl, nd, &iwork[inode - 1], &iwork[ndiml - 1], &iwork[ndimr - 1], smlsiz);
    //
    //     The following code applies back the left singular vector factors.
    //     For applying back the right singular vector factors, go to 50.
    //
    if (icompq == 1) {
        goto statement_50;
    }
    //
    //     The nodes on the bottom level of the tree were solved
    //     by Rlasdq. The corresponding left and right singular vector
    //     matrices are in explicit form. First apply back the left
    //     singular vector matrices.
    //
    ndb1 = (nd + 1) / 2;
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
        nr = iwork[(ndimr + i1) - 1];
        nlf = ic - nl;
        nrf = ic + 1;
        Rgemm("T", "N", nl, nrhs, nl, one, &u[(nlf - 1)], ldu, &b[(nlf - 1)], ldb, zero, &bx[(nlf - 1)], ldbx);
        Rgemm("T", "N", nr, nrhs, nr, one, &u[(nrf - 1)], ldu, &b[(nrf - 1)], ldb, zero, &bx[(nrf - 1)], ldbx);
    }
    //
    //     Next copy the rows of B that correspond to unchanged rows
    //     in the bidiagonal matrix to BX.
    //
    for (i = 1; i <= nd; i = i + 1) {
        ic = iwork[(inode + i - 1) - 1];
        Rcopy(nrhs, &b[(ic - 1)], ldb, &bx[(ic - 1)], ldbx);
    }
    //
    //     Finally go through the left singular vector matrices of all
    //     the other subproblems bottom-up on the tree.
    //
    j = pow(2, nlvl);
    sqre = 0;
    //
    for (lvl = nlvl; lvl >= 1; lvl = lvl - 1) {
        lvl2 = 2 * lvl - 1;
        //
        //        find the first node LF and last node LL on
        //        the current level LVL
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
            j = j - 1;
            Rlals0(icompq, nl, nr, sqre, nrhs, &bx[(nlf - 1)], ldbx, &b[(nlf - 1)], ldb, &perm[(nlf - 1) + (lvl - 1) * ldperm], givptr[j - 1], &givcol[(nlf - 1) + (lvl2 - 1) * ldgivcol], ldgcol, &givnum[(nlf - 1) + (lvl2 - 1) * ldgivnum], ldu, &poles[(nlf - 1) + (lvl2 - 1) * ldpoles], &difl[(nlf - 1) + (lvl - 1) * lddifl], &difr[(nlf - 1) + (lvl2 - 1) * lddifr], &z[(nlf - 1) + (lvl - 1) * ldz], k[j - 1], c[j - 1], s[j - 1], work, info);
        }
    }
    goto statement_90;
//
//     ICOMPQ = 1: applying back the right singular vector factors.
//
statement_50:
    //
    //     First now go through the right singular vector matrices of all
    //     the tree nodes top-down.
    //
    j = 0;
    for (lvl = 1; lvl <= nlvl; lvl = lvl + 1) {
        lvl2 = 2 * lvl - 1;
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
        for (i = ll; i >= lf; i = i - 1) {
            im1 = i - 1;
            ic = iwork[(inode + im1) - 1];
            nl = iwork[(ndiml + im1) - 1];
            nr = iwork[(ndimr + im1) - 1];
            nlf = ic - nl;
            nrf = ic + 1;
            if (i == ll) {
                sqre = 0;
            } else {
                sqre = 1;
            }
            j++;
            Rlals0(icompq, nl, nr, sqre, nrhs, &b[(nlf - 1)], ldb, &bx[(nlf - 1)], ldbx, &perm[(nlf - 1) + (lvl - 1) * ldperm], givptr[j - 1], &givcol[(nlf - 1) + (lvl2 - 1) * ldgivcol], ldgcol, &givnum[(nlf - 1) + (lvl2 - 1) * ldgivnum], ldu, &poles[(nlf - 1) + (lvl2 - 1) * ldpoles], &difl[(nlf - 1) + (lvl - 1) * lddifl], &difr[(nlf - 1) + (lvl2 - 1) * lddifr], &z[(nlf - 1) + (lvl - 1) * ldz], k[j - 1], c[j - 1], s[j - 1], work, info);
        }
    }
    //
    //     The nodes on the bottom level of the tree were solved
    //     by Rlasdq. The corresponding right singular vector
    //     matrices are in explicit form. Apply them back.
    //
    ndb1 = (nd + 1) / 2;
    for (i = ndb1; i <= nd; i = i + 1) {
        i1 = i - 1;
        ic = iwork[(inode + i1) - 1];
        nl = iwork[(ndiml + i1) - 1];
        nr = iwork[(ndimr + i1) - 1];
        nlp1 = nl + 1;
        if (i == nd) {
            nrp1 = nr;
        } else {
            nrp1 = nr + 1;
        }
        nlf = ic - nl;
        nrf = ic + 1;
        Rgemm("T", "N", nlp1, nrhs, nlp1, one, &vt[(nlf - 1)], ldu, &b[(nlf - 1)], ldb, zero, &bx[(nlf - 1)], ldbx);
        Rgemm("T", "N", nrp1, nrhs, nrp1, one, &vt[(nrf - 1)], ldu, &b[(nrf - 1)], ldb, zero, &bx[(nrf - 1)], ldbx);
    }
//
statement_90:;
    //
    //     End of Rlalsa
    //
}

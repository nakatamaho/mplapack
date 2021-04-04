/*
 * Copyright (c) 2021
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

void Clalsa(INTEGER const &icompq, INTEGER const &smlsiz, INTEGER const &n, INTEGER const &nrhs, COMPLEX *b, INTEGER const &ldb, COMPLEX *bx, INTEGER const &ldbx, REAL *u, INTEGER const &ldu, REAL *vt, INTEGER *k, REAL *difl, REAL *difr, REAL *z, REAL *poles, INTEGER *givptr, arr_cref<INTEGER, 2> givcol, INTEGER const &ldgcol, arr_cref<INTEGER, 2> perm, REAL *givnum, REAL *c, REAL *s, REAL *rwork, INTEGER *iwork, INTEGER &info) {
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
    INTEGER j = 0;
    INTEGER jcol = 0;
    INTEGER jrow = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    INTEGER jreal = 0;
    INTEGER jimag = 0;
    INTEGER sqre = 0;
    INTEGER lvl = 0;
    INTEGER lvl2 = 0;
    INTEGER lf = 0;
    INTEGER ll = 0;
    INTEGER im1 = 0;
    INTEGER nlp1 = 0;
    INTEGER nrp1 = 0;
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
    //     .. Intrinsic Functions ..
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
        Mxerbla("Clalsa", -info);
        return;
    }
    //
    //     Book-keeping and  setting up the computation tree.
    //
    inode = 1;
    ndiml = inode + n;
    ndimr = ndiml + n;
    //
    Rlasdt(n, nlvl, nd, iwork[inode - 1], iwork[ndiml - 1], iwork[ndimr - 1], smlsiz);
    //
    //     The following code applies back the left singular vector factors.
    //     For applying back the right singular vector factors, go to 170.
    //
    if (icompq == 1) {
        goto statement_170;
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
        //
        //        Since B and BX are complex, the following call to Rgemm
        //        is performed in two steps (real and imaginary parts).
        //
        //        CALL Rgemm( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
        //     $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
        //
        j = nl * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nlf; jrow <= nlf + nl - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].real();
            }
        }
        Rgemm("T", "N", nl, nrhs, nl, one, u[(nlf - 1)], ldu, rwork[(1 + nl * nrhs * 2) - 1], nl, zero, rwork[1 - 1], nl);
        j = nl * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nlf; jrow <= nlf + nl - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].imag();
            }
        }
        Rgemm("T", "N", nl, nrhs, nl, one, u[(nlf - 1)], ldu, rwork[(1 + nl * nrhs * 2) - 1], nl, zero, rwork[(1 + nl * nrhs) - 1], nl);
        jreal = 0;
        jimag = nl * nrhs;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nlf; jrow <= nlf + nl - 1; jrow = jrow + 1) {
                jreal++;
                jimag++;
                bx[(jrow - 1) + (jcol - 1) * ldbx] = COMPLEX(rwork[jreal - 1], rwork[jimag - 1]);
            }
        }
        //
        //        Since B and BX are complex, the following call to Rgemm
        //        is performed in two steps (real and imaginary parts).
        //
        //        CALL Rgemm( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
        //    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
        //
        j = nr * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nrf; jrow <= nrf + nr - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].real();
            }
        }
        Rgemm("T", "N", nr, nrhs, nr, one, u[(nrf - 1)], ldu, rwork[(1 + nr * nrhs * 2) - 1], nr, zero, rwork[1 - 1], nr);
        j = nr * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nrf; jrow <= nrf + nr - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].imag();
            }
        }
        Rgemm("T", "N", nr, nrhs, nr, one, u[(nrf - 1)], ldu, rwork[(1 + nr * nrhs * 2) - 1], nr, zero, rwork[(1 + nr * nrhs) - 1], nr);
        jreal = 0;
        jimag = nr * nrhs;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nrf; jrow <= nrf + nr - 1; jrow = jrow + 1) {
                jreal++;
                jimag++;
                bx[(jrow - 1) + (jcol - 1) * ldbx] = COMPLEX(rwork[jreal - 1], rwork[jimag - 1]);
            }
        }
        //
    }
    //
    //     Next copy the rows of B that correspond to unchanged rows
    //     in the bidiagonal matrix to BX.
    //
    for (i = 1; i <= nd; i = i + 1) {
        ic = iwork[(inode + i - 1) - 1];
        Ccopy(nrhs, b[(ic - 1)], ldb, bx[(ic - 1)], ldbx);
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
            Clals0(icompq, nl, nr, sqre, nrhs, bx[(nlf - 1)], ldbx, b[(nlf - 1)], ldb, perm[(nlf - 1) + (lvl - 1) * ldperm], givptr[j - 1], givcol[(nlf - 1) + (lvl2 - 1) * ldgivcol], ldgcol, givnum[(nlf - 1) + (lvl2 - 1) * ldgivnum], ldu, poles[(nlf - 1) + (lvl2 - 1) * ldpoles], difl[(nlf - 1) + (lvl - 1) * lddifl], difr[(nlf - 1) + (lvl2 - 1) * lddifr], z[(nlf - 1) + (lvl - 1) * ldz], k[j - 1], c[j - 1], s[j - 1], rwork, info);
        }
    }
    goto statement_330;
//
//     ICOMPQ = 1: applying back the right singular vector factors.
//
statement_170:
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
            Clals0(icompq, nl, nr, sqre, nrhs, b[(nlf - 1)], ldb, bx[(nlf - 1)], ldbx, perm[(nlf - 1) + (lvl - 1) * ldperm], givptr[j - 1], givcol[(nlf - 1) + (lvl2 - 1) * ldgivcol], ldgcol, givnum[(nlf - 1) + (lvl2 - 1) * ldgivnum], ldu, poles[(nlf - 1) + (lvl2 - 1) * ldpoles], difl[(nlf - 1) + (lvl - 1) * lddifl], difr[(nlf - 1) + (lvl2 - 1) * lddifr], z[(nlf - 1) + (lvl - 1) * ldz], k[j - 1], c[j - 1], s[j - 1], rwork, info);
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
        //
        //        Since B and BX are complex, the following call to Rgemm is
        //        performed in two steps (real and imaginary parts).
        //
        //        CALL Rgemm( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
        //    $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
        //
        j = nlp1 * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nlf; jrow <= nlf + nlp1 - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].real();
            }
        }
        Rgemm("T", "N", nlp1, nrhs, nlp1, one, vt[(nlf - 1)], ldu, rwork[(1 + nlp1 * nrhs * 2) - 1], nlp1, zero, rwork[1 - 1], nlp1);
        j = nlp1 * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nlf; jrow <= nlf + nlp1 - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].imag();
            }
        }
        Rgemm("T", "N", nlp1, nrhs, nlp1, one, vt[(nlf - 1)], ldu, rwork[(1 + nlp1 * nrhs * 2) - 1], nlp1, zero, rwork[(1 + nlp1 * nrhs) - 1], nlp1);
        jreal = 0;
        jimag = nlp1 * nrhs;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nlf; jrow <= nlf + nlp1 - 1; jrow = jrow + 1) {
                jreal++;
                jimag++;
                bx[(jrow - 1) + (jcol - 1) * ldbx] = COMPLEX(rwork[jreal - 1], rwork[jimag - 1]);
            }
        }
        //
        //        Since B and BX are complex, the following call to Rgemm is
        //        performed in two steps (real and imaginary parts).
        //
        //        CALL Rgemm( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
        //    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
        //
        j = nrp1 * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nrf; jrow <= nrf + nrp1 - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].real();
            }
        }
        Rgemm("T", "N", nrp1, nrhs, nrp1, one, vt[(nrf - 1)], ldu, rwork[(1 + nrp1 * nrhs * 2) - 1], nrp1, zero, rwork[1 - 1], nrp1);
        j = nrp1 * nrhs * 2;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nrf; jrow <= nrf + nrp1 - 1; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].imag();
            }
        }
        Rgemm("T", "N", nrp1, nrhs, nrp1, one, vt[(nrf - 1)], ldu, rwork[(1 + nrp1 * nrhs * 2) - 1], nrp1, zero, rwork[(1 + nrp1 * nrhs) - 1], nrp1);
        jreal = 0;
        jimag = nrp1 * nrhs;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = nrf; jrow <= nrf + nrp1 - 1; jrow = jrow + 1) {
                jreal++;
                jimag++;
                bx[(jrow - 1) + (jcol - 1) * ldbx] = COMPLEX(rwork[jreal - 1], rwork[jimag - 1]);
            }
        }
        //
    }
//
statement_330:;
    //
    //     End of Clalsa
    //
}

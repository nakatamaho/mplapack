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

void Rlals0(INTEGER const icompq, INTEGER const nl, INTEGER const nr, INTEGER const sqre, INTEGER const nrhs, REAL *b, INTEGER const ldb, REAL *bx, INTEGER const ldbx, INTEGER *perm, INTEGER const givptr, INTEGER *givcol, INTEGER const ldgcol, REAL *givnum, INTEGER const ldgnum, REAL *poles, REAL *difl, REAL *difr, REAL *z, INTEGER const k, REAL const c, REAL const s, REAL *work, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    info = 0;
    INTEGER n = nl + nr + 1;
    //
    if ((icompq < 0) || (icompq > 1)) {
        info = -1;
    } else if (nl < 1) {
        info = -2;
    } else if (nr < 1) {
        info = -3;
    } else if ((sqre < 0) || (sqre > 1)) {
        info = -4;
    } else if (nrhs < 1) {
        info = -5;
    } else if (ldb < n) {
        info = -7;
    } else if (ldbx < n) {
        info = -9;
    } else if (givptr < 0) {
        info = -11;
    } else if (ldgcol < n) {
        info = -13;
    } else if (ldgnum < n) {
        info = -15;
    } else if (k < 1) {
        info = -20;
    }
    if (info != 0) {
        Mxerbla("Rlals0", -info);
        return;
    }
    //
    INTEGER m = n + sqre;
    INTEGER nlp1 = nl + 1;
    //
    INTEGER i = 0;
    const REAL zero = 0.0;
    const REAL negone = -1.0;
    INTEGER j = 0;
    REAL diflj = 0.0;
    REAL dj = 0.0;
    REAL dsigj = 0.0;
    REAL difrj = 0.0;
    REAL dsigjp = 0.0;
    REAL temp = 0.0;
    const REAL one = 1.0;
    if (icompq == 0) {
        //
        //        Apply back orthogonal transformations from the left.
        //
        //        Step (1L): apply back the Givens rotations performed.
        //
        for (i = 1; i <= givptr; i = i + 1) {
            Rrot(nrhs, &b[(givcol[(i - 1) + (2 - 1) * ldgcol] - 1)], ldb, &b[(givcol[(i - 1)] - 1)], ldb, givnum[(i - 1) + (2 - 1) * ldgnum], givnum[(i - 1)]);
        }
        //
        //        Step (2L): permute rows of B.
        //
        Rcopy(nrhs, &b[(nlp1 - 1)], ldb, &bx[(1 - 1)], ldbx);
        for (i = 2; i <= n; i = i + 1) {
            Rcopy(nrhs, &b[(perm[i - 1] - 1)], ldb, &bx[(i - 1)], ldbx);
        }
        //
        //        Step (3L): apply the inverse of the left singular vector
        //        matrix to BX.
        //
        if (k == 1) {
            Rcopy(nrhs, bx, ldbx, b, ldb);
            if (z[1 - 1] < zero) {
                Rscal(nrhs, negone, b, ldb);
            }
        } else {
            for (j = 1; j <= k; j = j + 1) {
                diflj = difl[j - 1];
                dj = poles[(j - 1)];
                dsigj = -poles[(j - 1) + (2 - 1) * ldgnum];
                if (j < k) {
                    difrj = -difr[(j - 1)];
                    dsigjp = -poles[((j + 1) - 1) + (2 - 1) * ldgnum];
                }
                if ((z[j - 1] == zero) || (poles[(j - 1) + (2 - 1) * ldgnum] == zero)) {
                    work[j - 1] = zero;
                } else {
                    work[j - 1] = -poles[(j - 1) + (2 - 1) * ldgnum] * z[j - 1] / diflj / (poles[(j - 1) + (2 - 1) * ldgnum] + dj);
                }
                for (i = 1; i <= j - 1; i = i + 1) {
                    if ((z[i - 1] == zero) || (poles[(i - 1) + (2 - 1) * ldgnum] == zero)) {
                        work[i - 1] = zero;
                    } else {
                        work[i - 1] = poles[(i - 1) + (2 - 1) * ldgnum] * z[i - 1] / (Rlamc3(poles[(i - 1) + (2 - 1) * ldgnum], dsigj) - diflj) / (poles[(i - 1) + (2 - 1) * ldgnum] + dj);
                    }
                }
                for (i = j + 1; i <= k; i = i + 1) {
                    if ((z[i - 1] == zero) || (poles[(i - 1) + (2 - 1) * ldgnum] == zero)) {
                        work[i - 1] = zero;
                    } else {
                        work[i - 1] = poles[(i - 1) + (2 - 1) * ldgnum] * z[i - 1] / (Rlamc3(poles[(i - 1) + (2 - 1) * ldgnum], dsigjp) + difrj) / (poles[(i - 1) + (2 - 1) * ldgnum] + dj);
                    }
                }
                work[1 - 1] = negone;
                temp = Rnrm2(k, work, 1);
                Rgemv("T", k, nrhs, one, bx, ldbx, work, 1, zero, &b[(j - 1)], ldb);
                Rlascl("G", 0, 0, temp, one, 1, nrhs, &b[(j - 1)], ldb, info);
            }
        }
        //
        //        Move the deflated rows of BX to B also.
        //
        if (k < max(m, n)) {
            Rlacpy("A", n - k, nrhs, &bx[((k + 1) - 1)], ldbx, &b[((k + 1) - 1)], ldb);
        }
    } else {
        //
        //        Apply back the right orthogonal transformations.
        //
        //        Step (1R): apply back the new right singular vector matrix
        //        to B.
        //
        if (k == 1) {
            Rcopy(nrhs, b, ldb, bx, ldbx);
        } else {
            for (j = 1; j <= k; j = j + 1) {
                dsigj = poles[(j - 1) + (2 - 1) * ldgnum];
                if (z[j - 1] == zero) {
                    work[j - 1] = zero;
                } else {
                    work[j - 1] = -z[j - 1] / difl[j - 1] / (dsigj + poles[(j - 1)]) / difr[(j - 1) + (2 - 1) * ldgnum];
                }
                for (i = 1; i <= j - 1; i = i + 1) {
                    if (z[j - 1] == zero) {
                        work[i - 1] = zero;
                    } else {
                        work[i - 1] = z[j - 1] / (Rlamc3(dsigj, -poles[((i + 1) - 1) + (2 - 1) * ldgnum]) - difr[(i - 1)]) / (dsigj + poles[(i - 1)]) / difr[(i - 1) + (2 - 1) * ldgnum];
                    }
                }
                for (i = j + 1; i <= k; i = i + 1) {
                    if (z[j - 1] == zero) {
                        work[i - 1] = zero;
                    } else {
                        work[i - 1] = z[j - 1] / (Rlamc3(dsigj, -poles[(i - 1) + (2 - 1) * ldgnum]) - difl[i - 1]) / (dsigj + poles[(i - 1)]) / difr[(i - 1) + (2 - 1) * ldgnum];
                    }
                }
                Rgemv("T", k, nrhs, one, b, ldb, work, 1, zero, &bx[(j - 1)], ldbx);
            }
        }
        //
        //        Step (2R): if SQRE = 1, apply back the rotation that is
        //        related to the right null space of the subproblem.
        //
        if (sqre == 1) {
            Rcopy(nrhs, &b[(m - 1)], ldb, &bx[(m - 1)], ldbx);
            Rrot(nrhs, &bx[(1 - 1)], ldbx, &bx[(m - 1)], ldbx, c, s);
        }
        if (k < max(m, n)) {
            Rlacpy("A", n - k, nrhs, &b[((k + 1) - 1)], ldb, &bx[((k + 1) - 1)], ldbx);
        }
        //
        //        Step (3R): permute rows of B.
        //
        Rcopy(nrhs, &bx[(1 - 1)], ldbx, &b[(nlp1 - 1)], ldb);
        if (sqre == 1) {
            Rcopy(nrhs, &bx[(m - 1)], ldbx, &b[(m - 1)], ldb);
        }
        for (i = 2; i <= n; i = i + 1) {
            Rcopy(nrhs, &bx[(i - 1)], ldbx, &b[(perm[i - 1] - 1)], ldb);
        }
        //
        //        Step (4R): apply back the Givens rotations performed.
        //
        for (i = givptr; i >= 1; i = i - 1) {
            Rrot(nrhs, &b[(givcol[(i - 1) + (2 - 1) * ldgcol] - 1)], ldb, &b[(givcol[(i - 1)] - 1)], ldb, givnum[(i - 1) + (2 - 1) * ldgnum], -givnum[(i - 1)]);
        }
    }
    //
    //     End of Rlals0
    //
}

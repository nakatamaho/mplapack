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

// Derived from LAPACK routine ZLATRS3.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Clatrs3(const char *uplo, const char *trans, const char *diag, const char *normin, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const ldx, REAL *scale, REAL *cnorm, REAL *work, INTEGER const lwork, INTEGER &info) {
    //
    // .. Scalar Arguments ..
    // ..
    // .. Array Arguments ..
    // ..
    //
    // =====================================================================
    //
    // .. Parameters ..
    // ..
    // .. Local Arrays ..
    // ..
    // .. Local Scalars ..
    // ..
    // .. External Functions ..
    // ..
    // .. External Subroutines ..
    // ..
    // .. Intrinsic Functions ..
    // ..
    // .. Executable Statements ..
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    bool notran = Mlsame(trans, "N");
    bool nounit = Mlsame(diag, "N");
    bool lquery = (lwork == -1);
    //
    // Partition A and X into blocks.
    //
    const INTEGER nbmin = 8;
    INTEGER nb = max(nbmin, iMlaenv(1, "Clatrs", "", n, n, -1, -1));
    const INTEGER nbmax = 64;
    nb = min(nbmax, nb);
    INTEGER nba = max((INTEGER)1, (n + nb - 1) / nb);
    const INTEGER nbrhs = 32;
    INTEGER nbx = max((INTEGER)1, (nrhs + nbrhs - 1) / nbrhs);
    //
    // Compute the workspace
    //
    // The workspace comprises two parts.
    // The first part stores the local scale factors. Each simultaneously
    // computed right-hand side requires one local scale factor per block
    // row. WORK( I + KK * LDS ) is the scale factor of the vector
    // segment associated with the I-th block row and the KK-th vector
    // in the block column.
    //
    INTEGER lscale = nba * max(nba, min(nrhs, nbrhs));
    INTEGER lds = nba;
    //
    // The second part stores upper bounds of the triangular A. There are
    // a total of NBA x NBA blocks, of which only the upper triangular
    // part or the lower triangular part is referenced. The upper bound of
    // the block A( I, J ) is stored as WORK( AWRK + I + J * NBA ).
    //
    INTEGER lanrm = nba * nba;
    INTEGER awrk = lscale;
    //
    INTEGER lwmin = 0;
    if (min(n, nrhs) == 0) {
        lwmin = 1;
    } else {
        lwmin = lscale + lanrm;
    }
    work[1 - 1] = lwmin;
    //
    // Test the input parameters.
    //
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = -2;
    } else if (!nounit && !Mlsame(diag, "U")) {
        info = -3;
    } else if (!Mlsame(normin, "Y") && !Mlsame(normin, "N")) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (nrhs < 0) {
        info = -6;
    } else if (lda < max((INTEGER)1, n)) {
        info = -8;
    } else if (ldx < max((INTEGER)1, n)) {
        info = -10;
    } else if (!lquery && lwork < lwmin) {
        info = -14;
    }
    if (info != 0) {
        Mxerbla("Clatrs3", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    // Initialize scaling factors
    //
    INTEGER kk = 0;
    const REAL one = 1.0;
    for (kk = 1; kk <= nrhs; kk = kk + 1) {
        scale[kk - 1] = one;
    }
    //
    // Quick return if possible
    //
    if (min(n, nrhs) == 0) {
        return;
    }
    //
    // Determine machine dependent constant to control overflow.
    //
    REAL bignum = Rlamch("Overflow");
    REAL smlnum = Rlamch("Safe Minimum");
    //
    // Use unblocked code for small problems
    //
    const INTEGER nrhsmin = 2;
    INTEGER k = 0;
    if (nrhs < nrhsmin) {
        Clatrs(uplo, trans, diag, normin, n, a, lda, &x[0], scale[1 - 1], cnorm, info);
        for (k = 2; k <= nrhs; k = k + 1) {
            Clatrs(uplo, trans, diag, "Y", n, a, lda, &x[(k - 1) * ldx], scale[k - 1], cnorm, info);
        }
        return;
    }
    //
    // Compute norms of blocks of A excluding diagonal blocks and find
    // the block with the largest norm TMAX.
    //
    const REAL zero = 0.0;
    REAL tmax = zero;
    INTEGER j = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    INTEGER ifirst = 0;
    INTEGER ilast = 0;
    INTEGER i = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    REAL w[nbmax];
    REAL anrm = 0.0;
    for (j = 1; j <= nba; j = j + 1) {
        j1 = (j - 1) * nb + 1;
        j2 = min(j * nb, n) + 1;
        if (upper) {
            ifirst = 1;
            ilast = j - 1;
        } else {
            ifirst = j + 1;
            ilast = nba;
        }
        for (i = ifirst; i <= ilast; i = i + 1) {
            i1 = (i - 1) * nb + 1;
            i2 = min(i * nb, n) + 1;
            //
            // Compute upper bound of A( I1:I2-1, J1:J2-1 ).
            //
            if (notran) {
                anrm = Clange("I", i2 - i1, j2 - j1, &a[(i1 - 1) + (j1 - 1) * lda], lda, w);
                work[(awrk + i + (j - 1) * nba) - 1] = anrm;
            } else {
                anrm = Clange("1", i2 - i1, j2 - j1, &a[(i1 - 1) + (j1 - 1) * lda], lda, w);
                work[(awrk + j + (i - 1) * nba) - 1] = anrm;
            }
            tmax = max(tmax, anrm);
        }
    }
    //
    if (!(tmax <= Rlamch("Overflow"))) {
        //
        // Some matrix entries have huge absolute value. At least one upper
        // bound norm( A(I1:I2-1, J1:J2-1), 'I') is not a valid floating-point
        // number, either due to overflow in LANGE or due to Inf in A.
        // Fall back to LATRS. Set normin = 'N' for every right-hand side to
        // force computation of TSCAL in LATRS to avoid the likely overflow
        // in the computation of the column norms CNORM.
        //
        for (k = 1; k <= nrhs; k = k + 1) {
            Clatrs(uplo, trans, diag, "N", n, a, lda, &x[(k - 1) * ldx], scale[k - 1], cnorm, info);
        }
        return;
    }
    //
    // Every right-hand side requires workspace to store NBA local scale
    // factors. To save workspace, X is computed successively in block columns
    // of width NBRHS, requiring a total of NBA x NBRHS space. If sufficient
    // workspace is available, larger values of NBRHS or NBRHS = NRHS are viable.
    INTEGER k1 = 0;
    INTEGER k2 = 0;
    INTEGER jfirst = 0;
    INTEGER jlast = 0;
    INTEGER jinc = 0;
    INTEGER rhs = 0;
    REAL scaloc = 0.0;
    REAL xnrm[nbrhs];
    INTEGER ii = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    REAL scal = 0.0;
    REAL rscal = 0.0;
    INTEGER iinc = 0;
    REAL scamin = 0.0;
    REAL bnrm = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    for (k = 1; k <= nbx; k = k + 1) {
        // Loop over block columns (index = K) of X and, for column-wise scalings,
        // over individual columns (index = KK).
        // K1: column index of the first column in X( J, K )
        // K2: column index of the first column in X( J, K+1 )
        // so the K2 - K1 is the column count of the block X( J, K )
        k1 = (k - 1) * nbrhs + 1;
        k2 = min(k * nbrhs, nrhs) + 1;
        //
        // Initialize local scaling factors of current block column X( J, K )
        //
        for (kk = 1; kk <= k2 - k1; kk = kk + 1) {
            for (i = 1; i <= nba; i = i + 1) {
                work[(i + kk * lds) - 1] = one;
            }
        }
        //
        if (notran) {
            //
            // Solve A * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))
            //
            if (upper) {
                jfirst = nba;
                jlast = 1;
                jinc = -1;
            } else {
                jfirst = 1;
                jlast = nba;
                jinc = 1;
            }
        } else {
            //
            // Solve op(A) * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))
            // where op(A) = A**T or op(A) = A**H
            //
            if (upper) {
                jfirst = 1;
                jlast = nba;
                jinc = 1;
            } else {
                jfirst = nba;
                jlast = 1;
                jinc = -1;
            }
        }
        //
        for (j = jfirst; jinc > 0 ? j <= jlast : j >= jlast; j = j + jinc) {
            // J1: row index of the first row in A( J, J )
            // J2: row index of the first row in A( J+1, J+1 )
            // so that J2 - J1 is the row count of the block A( J, J )
            j1 = (j - 1) * nb + 1;
            j2 = min(j * nb, n) + 1;
            //
            // Solve op(A( J, J )) * X( J, RHS ) = SCALOC * B( J, RHS )
            //
            for (kk = 1; kk <= k2 - k1; kk = kk + 1) {
                rhs = k1 + kk - 1;
                if (kk == 1) {
                    Clatrs(uplo, trans, diag, "N", j2 - j1, &a[(j1 - 1) + (j1 - 1) * lda], lda, &x[(j1 - 1) + (rhs - 1) * ldx], scaloc, cnorm, info);
                } else {
                    Clatrs(uplo, trans, diag, "Y", j2 - j1, &a[(j1 - 1) + (j1 - 1) * lda], lda, &x[(j1 - 1) + (rhs - 1) * ldx], scaloc, cnorm, info);
                }
                // Find largest absolute value entry in the vector segment
                // X( J1:J2-1, RHS ) as an upper bound for the worst case
                // growth in the linear updates.
                xnrm[kk - 1] = Clange("I", j2 - j1, 1, &x[(j1 - 1) + (rhs - 1) * ldx], ldx, w);
                //
                if (scaloc == zero) {
                    // LATRS found that A is singular through A(j,j) = 0.
                    // Reset the computation x(1:n) = 0, x(j) = 1, SCALE = 0
                    // and compute op(A)*x = 0. Note that X(J1:J2-1, KK) is
                    // set by LATRS.
                    scale[rhs - 1] = zero;
                    for (ii = 1; ii <= j1 - 1; ii = ii + 1) {
                        x[(ii - 1) + (kk - 1) * ldx] = czero;
                    }
                    for (ii = j2; ii <= n; ii = ii + 1) {
                        x[(ii - 1) + (kk - 1) * ldx] = czero;
                    }
                    // Discard the local scale factors.
                    for (ii = 1; ii <= nba; ii = ii + 1) {
                        work[(ii + kk * lds) - 1] = one;
                    }
                    scaloc = one;
                } else if (scaloc * work[(j + kk * lds) - 1] == zero) {
                    // LATRS computed a valid scale factor, but combined with
                    // the current scaling the solution does not have a
                    // scale factor > 0.
                    //
                    // Set WORK( J+KK*LDS ) to smallest valid scale
                    // factor and increase SCALOC accordingly.
                    scal = work[(j + kk * lds) - 1] / smlnum;
                    scaloc = scaloc * scal;
                    work[(j + kk * lds) - 1] = smlnum;
                    // If LATRS overestimated the growth, x may be
                    // rescaled to preserve a valid combined scale
                    // factor WORK( J, KK ) > 0.
                    rscal = one / scaloc;
                    if (xnrm[kk - 1] * rscal <= bignum) {
                        xnrm[kk - 1] = xnrm[kk - 1] * rscal;
                        CRscal(j2 - j1, rscal, &x[(j1 - 1) + (rhs - 1) * ldx], 1);
                        scaloc = one;
                    } else {
                        // The system op(A) * x = b is badly scaled and its
                        // solution cannot be represented as (1/scale) * x.
                        // Set x to zero. This approach deviates from LATRS
                        // where a completely meaningless non-zero vector
                        // is returned that is not a solution to op(A) * x = b.
                        scale[rhs - 1] = zero;
                        for (ii = 1; ii <= n; ii = ii + 1) {
                            x[(ii - 1) + (kk - 1) * ldx] = czero;
                        }
                        // Discard the local scale factors.
                        for (ii = 1; ii <= nba; ii = ii + 1) {
                            work[(ii + kk * lds) - 1] = one;
                        }
                        scaloc = one;
                    }
                }
                scaloc = scaloc * work[(j + kk * lds) - 1];
                work[(j + kk * lds) - 1] = scaloc;
            }
            //
            // Linear block updates
            //
            if (notran) {
                if (upper) {
                    ifirst = j - 1;
                    ilast = 1;
                    iinc = -1;
                } else {
                    ifirst = j + 1;
                    ilast = nba;
                    iinc = 1;
                }
            } else {
                if (upper) {
                    ifirst = j + 1;
                    ilast = nba;
                    iinc = 1;
                } else {
                    ifirst = j - 1;
                    ilast = 1;
                    iinc = -1;
                }
            }
            //
            for (i = ifirst; iinc > 0 ? i <= ilast : i >= ilast; i = i + iinc) {
                // I1: row index of the first column in X( I, K )
                // I2: row index of the first column in X( I+1, K )
                // so the I2 - I1 is the row count of the block X( I, K )
                i1 = (i - 1) * nb + 1;
                i2 = min(i * nb, n) + 1;
                //
                // Prepare the linear update to be executed with GEMM.
                // For each column, compute a consistent scaling, a
                // scaling factor to survive the linear update, and
                // rescale the column segments, if necessary. Then
                // the linear update is safely executed.
                //
                for (kk = 1; kk <= k2 - k1; kk = kk + 1) {
                    rhs = k1 + kk - 1;
                    // Compute consistent scaling
                    scamin = min(work[(i + kk * lds) - 1], work[(j + kk * lds) - 1]);
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    bnrm = Clange("I", i2 - i1, 1, &x[(i1 - 1) + (rhs - 1) * ldx], ldx, w);
                    bnrm = bnrm * (scamin / work[(i + kk * lds) - 1]);
                    xnrm[kk - 1] = xnrm[kk - 1] * (scamin / work[(j + kk * lds) - 1]);
                    anrm = work[(awrk + i + (j - 1) * nba) - 1];
                    scaloc = Rlarmm(anrm, xnrm[kk - 1], bnrm);
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to X( I, KK ) and X( J, KK ).
                    //
                    scal = (scamin / work[(i + kk * lds) - 1]) * scaloc;
                    if (scal != one) {
                        CRscal(i2 - i1, scal, &x[(i1 - 1) + (rhs - 1) * ldx], 1);
                        work[(i + kk * lds) - 1] = scamin * scaloc;
                    }
                    //
                    scal = (scamin / work[(j + kk * lds) - 1]) * scaloc;
                    if (scal != one) {
                        CRscal(j2 - j1, scal, &x[(j1 - 1) + (rhs - 1) * ldx], 1);
                        work[(j + kk * lds) - 1] = scamin * scaloc;
                    }
                }
                //
                if (notran) {
                    //
                    // B( I, K ) := B( I, K ) - A( I, J ) * X( J, K )
                    //
                    Cgemm("N", "N", i2 - i1, k2 - k1, j2 - j1, -cone, &a[(i1 - 1) + (j1 - 1) * lda], lda, &x[(j1 - 1) + (k1 - 1) * ldx], ldx, cone, &x[(i1 - 1) + (k1 - 1) * ldx], ldx);
                } else if (Mlsame(trans, "T")) {
                    //
                    // B( I, K ) := B( I, K ) - A( I, J )**T * X( J, K )
                    //
                    Cgemm("T", "N", i2 - i1, k2 - k1, j2 - j1, -cone, &a[(j1 - 1) + (i1 - 1) * lda], lda, &x[(j1 - 1) + (k1 - 1) * ldx], ldx, cone, &x[(i1 - 1) + (k1 - 1) * ldx], ldx);
                } else {
                    //
                    // B( I, K ) := B( I, K ) - A( I, J )**H * X( J, K )
                    //
                    Cgemm("C", "N", i2 - i1, k2 - k1, j2 - j1, -cone, &a[(j1 - 1) + (i1 - 1) * lda], lda, &x[(j1 - 1) + (k1 - 1) * ldx], ldx, cone, &x[(i1 - 1) + (k1 - 1) * ldx], ldx);
                }
            }
        }
        //
        // Reduce local scaling factors
        //
        for (kk = 1; kk <= k2 - k1; kk = kk + 1) {
            rhs = k1 + kk - 1;
            for (i = 1; i <= nba; i = i + 1) {
                scale[rhs - 1] = min(scale[rhs - 1], work[(i + kk * lds) - 1]);
            }
        }
        //
        // Realize consistent scaling
        //
        for (kk = 1; kk <= k2 - k1; kk = kk + 1) {
            rhs = k1 + kk - 1;
            if (scale[rhs - 1] != one && scale[rhs - 1] != zero) {
                for (i = 1; i <= nba; i = i + 1) {
                    i1 = (i - 1) * nb + 1;
                    i2 = min(i * nb, n) + 1;
                    scal = scale[rhs - 1] / work[(i + kk * lds) - 1];
                    if (scal != one) {
                        CRscal(i2 - i1, scal, &x[(i1 - 1) + (rhs - 1) * ldx], 1);
                    }
                }
            }
        }
    }
    //
    // End of Clatrs3
    //
}

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

// Derived from LAPACK routine ZGEQP3RK.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Cgeqp3rk(INTEGER const m, INTEGER const n, INTEGER const nrhs, INTEGER const kmax, REAL &abstol, REAL &reltol, COMPLEX *a, INTEGER const lda, INTEGER &k, REAL &maxc2nrmk, REAL &relmaxc2nrmk, INTEGER *jpiv, COMPLEX *tau, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER &info) {
    //
    // Test input arguments
    // ====================
    //
    info = 0;
    bool lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (kmax < 0) {
        info = -4;
    } else if (Risnan(abstol)) {
        info = -5;
    } else if (Risnan(reltol)) {
        info = -6;
    } else if (lda < max((INTEGER)1, m)) {
        info = -8;
    }
    //
    // If the input parameters M, N, NRHS, KMAX, LDA are valid:
    // a) Test the input workspace size LWORK for the minimum
    // size requirement IWS.
    // b) Determine the optimal block size NB and optimal
    // workspace size LWKOPT to be returned in WORK(1)
    // in case of (1) LWORK < IWS, (2) LQUERY = .TRUE.,
    // (3) when routine exits.
    // Here, IWS is the miminum workspace required for unblocked
    // code.
    //
    INTEGER minmn = 0;
    INTEGER iws = 0;
    INTEGER lwkopt = 0;
    const INTEGER inb = 1;
    INTEGER nb = 0;
    if (info == 0) {
        minmn = min(m, n);
        if (minmn == 0) {
            iws = 1;
            lwkopt = 1;
        } else {
            //
            // Minimal workspace size in case of using only unblocked
            // BLAS 2 code in Claqp2rk.
            // 1) Claqp2rk: N+NRHS-1 to use in WORK array that is used
            // in Clarf1f subroutine inside Claqp2rk to apply an
            // elementary reflector from the left.
            // TOTAL_WORK_SIZE = 3*N + NRHS - 1
            //
            iws = n + nrhs - 1;
            //
            // Assign to NB optimal block size.
            //
            nb = iMlaenv(inb, "Cgeqp3rk", " ", m, n, -1, -1);
            //
            // A formula for the optimal workspace size in case of using
            // both unblocked BLAS 2 in Claqp2rk and blocked BLAS 3 code
            // in Claqp3rk.
            // 1) Cgeqp3rk, Claqp2rk, Claqp3rk: 2*N to store full and
            // partial column 2-norms.
            // 2) Claqp2rk: N+NRHS-1 to use in WORK array that is used
            // in Clarf1f subroutine to apply an elementary reflector
            // from the left.
            // 3) Claqp3rk: NB*(N+NRHS) to use in the work array F that
            // is used to apply a block reflector from
            // the left.
            // 4) Claqp3rk: NB to use in the auxilixary array AUX.
            // Sizes (2) and ((3) + (4)) should intersect, therefore
            // TOTAL_WORK_SIZE = 2*N + NB*( N+NRHS+1 ), given NBMIN=2.
            //
            lwkopt = 2 * n + nb * (n + nrhs + 1);
        }
        work[1 - 1] = COMPLEX(lwkopt);
        //
        if ((lwork < iws) && !lquery) {
            info = -15;
        }
    }
    //
    // NOTE: The optimal workspace size is returned in WORK(1), if
    // the input parameters M, N, NRHS, KMAX, LDA are valid.
    //
    if (info != 0) {
        Mxerbla("Cgeqp3rk", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    // Quick return if possible for M=0 or N=0.
    //
    const REAL zero = 0.0;
    if (minmn == 0) {
        k = 0;
        maxc2nrmk = zero;
        relmaxc2nrmk = zero;
        work[1 - 1] = COMPLEX(lwkopt);
        return;
    }
    //
    // ==================================================================
    //
    // Initialize column pivot array JPIV.
    //
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        jpiv[j - 1] = j;
    }
    //
    // ==================================================================
    //
    // Initialize storage for partial and exact column 2-norms.
    // a) The elements WORK(1:N) are used to store partial column
    // 2-norms of the matrix A, and may decrease in each computation
    // step; initialize to the values of complete columns 2-norms.
    // b) The elements WORK(N+1:2*N) are used to store complete column
    // 2-norms of the matrix A, they are not changed during the
    // computation; initialize the values of complete columns 2-norms.
    //
    for (j = 1; j <= n; j = j + 1) {
        rwork[j - 1] = RCnrm2(m, &a[(j - 1) * lda], 1);
        rwork[(n + j) - 1] = rwork[j - 1];
    }
    //
    // ==================================================================
    //
    // Compute the pivot column index and the maximum column 2-norm
    // for the whole original matrix stored in A(1:M,1:N).
    //
    INTEGER kp1 = iRamax(n, &rwork[1 - 1], 1);
    //
    // ==================================================================.
    //
    REAL maxc2nrm = 0.0;
    if (Risnan(maxc2nrm)) {
        //
        // Check if the matrix A contains NaN, set INFO parameter
        // to the column number where the first NaN is found and return
        // from the routine.
        //
        k = 0;
        info = kp1;
        //
        // Set MAXC2NRMK and  RELMAXC2NRMK to NaN.
        //
        maxc2nrmk = maxc2nrm;
        relmaxc2nrmk = maxc2nrm;
        //
        // Array TAU is not set and contains undefined elements.
        //
        work[1 - 1] = COMPLEX(lwkopt);
        return;
    }
    //
    // ===================================================================
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (maxc2nrm == zero) {
        //
        // Check is the matrix A is a zero matrix, set array TAU and
        // return from the routine.
        //
        k = 0;
        maxc2nrmk = zero;
        relmaxc2nrmk = zero;
        //
        for (j = 1; j <= minmn; j = j + 1) {
            tau[j - 1] = czero;
        }
        //
        work[1 - 1] = COMPLEX(lwkopt);
        return;
        //
    }
    //
    // ===================================================================
    //
    REAL hugeval = Rlamch("Overflow");
    //
    if (maxc2nrm > hugeval) {
        //
        // Check if the matrix A contains +Inf or -Inf, set INFO parameter
        // to the column number, where the first +/-Inf  is found plus N,
        // and continue the computation.
        //
        info = n + kp1;
        //
    }
    //
    // ==================================================================
    //
    // Quick return if possible for the case when the first
    // stopping criterion is satisfied, i.e. KMAX = 0.
    //
    const REAL one = 1.0;
    if (kmax == 0) {
        k = 0;
        maxc2nrmk = maxc2nrm;
        relmaxc2nrmk = one;
        for (j = 1; j <= minmn; j = j + 1) {
            tau[j - 1] = czero;
        }
        work[1 - 1] = COMPLEX(lwkopt);
        return;
    }
    //
    // ==================================================================
    //
    REAL eps = Rlamch("Epsilon");
    //
    // Adjust ABSTOL
    //
    REAL safmin = 0.0;
    const REAL two = 2.0;
    if (abstol >= zero) {
        safmin = Rlamch("Safe minimum");
        abstol = max(abstol, two * safmin);
    }
    //
    // Adjust RELTOL
    //
    if (reltol >= zero) {
        reltol = max(reltol, eps);
    }
    //
    // ===================================================================
    //
    // JMAX is the maximum index of the column to be factorized,
    // which is also limited by the first stopping criterion KMAX.
    //
    INTEGER jmax = min(kmax, minmn);
    //
    // ===================================================================
    //
    // Quick return if possible for the case when the second or third
    // stopping criterion for the whole original matrix is satified,
    // i.e. MAXC2NRM <= ABSTOL or RELMAXC2NRM <= RELTOL
    // (which is ONE <= RELTOL).
    //
    if (maxc2nrm <= abstol || one <= reltol) {
        //
        k = 0;
        maxc2nrmk = maxc2nrm;
        relmaxc2nrmk = one;
        //
        for (j = 1; j <= minmn; j = j + 1) {
            tau[j - 1] = czero;
        }
        //
        work[1 - 1] = COMPLEX(lwkopt);
        return;
    }
    //
    // ==================================================================
    // Factorize columns
    // ==================================================================
    //
    // Determine the block size.
    //
    INTEGER nbmin = 2;
    INTEGER nx = 0;
    //
    const INTEGER ixover = 3;
    const INTEGER inbmin = 2;
    if ((nb > 1) && (nb < minmn)) {
        //
        // Determine when to cross over from blocked to unblocked code.
        // (for N less than NX, unblocked code should be used).
        //
        nx = max((INTEGER)0, iMlaenv(ixover, "Cgeqp3rk", " ", m, n, -1, -1));
        //
        if (nx < minmn) {
            //
            // Determine if workspace is large enough for blocked code.
            //
            if (lwork < lwkopt) {
                //
                // Not enough workspace to use optimal block size that
                // is currently stored in NB.
                // Reduce NB and determine the minimum value of NB.
                //
                nb = (lwork - 2 * n) / (n + 1);
                nbmin = max((INTEGER)2, iMlaenv(inbmin, "Cgeqp3rk", " ", m, n, -1, -1));
                //
            }
        }
    }
    //
    // ==================================================================
    //
    // DONE is the boolean flag to rerpresent the case when the
    // factorization completed in the block factorization routine,
    // before the end of the block.
    //
    bool done = false;
    //
    // J is the column index.
    //
    j = 1;
    //
    // (1) Use blocked code initially.
    //
    // JMAXB is the maximum column index of the block, when the
    // blocked code is used, is also limited by the first stopping
    // criterion KMAX.
    //
    INTEGER jmaxb = min(kmax, minmn - nx);
    //
    INTEGER jb = 0;
    INTEGER n_sub = 0;
    INTEGER ioffset = 0;
    INTEGER jbf = 0;
    INTEGER iinfo = 0;
    if (nb >= nbmin && nb < jmax && jmaxb > 0) {
        //
        // Loop over the column blocks of the matrix A(1:M,1:JMAXB). Here:
        // J   is the column index of a column block;
        // JB  is the column block size to pass to block factorization
        // routine in a loop step;
        // JBF is the number of columns that were actually factorized
        // that was returned by the block factorization routine
        // in a loop step, JBF <= JB;
        // N_SUB is the number of columns in the submatrix;
        // IOFFSET is the number of rows that should not be factorized.
        //
        while (j <= jmaxb) {
            //
            jb = min(nb, jmaxb - j + 1);
            n_sub = n - j + 1;
            ioffset = j - 1;
            //
            // Factorize JB columns among the columns A(J:N).
            //
            Claqp3rk(m, n_sub, nrhs, ioffset, jb, abstol, reltol, kp1, maxc2nrm, &a[(j - 1) * lda], lda, done, jbf, maxc2nrmk, relmaxc2nrmk, &jpiv[j - 1], &tau[j - 1], &rwork[j - 1], &rwork[(n + j) - 1], &work[1 - 1], &work[(jb + 1) - 1], n + nrhs - j + 1, iwork, iinfo);
            //
            // Set INFO on the first occurence of Inf.
            //
            if (iinfo > n_sub && info == 0) {
                info = 2 * ioffset + iinfo;
            }
            //
            if (done) {
                //
                // Either the submatrix is zero before the end of the
                // column block, or ABSTOL or RELTOL criterion is
                // satisfied before the end of the column block, we can
                // return from the routine. Perform the following before
                // returning:
                // a) Set the number of factorized columns K,
                // K = IOFFSET + JBF from the last call of blocked
                // routine.
                // NOTE: 1) MAXC2NRMK and RELMAXC2NRMK are returned
                // by the block factorization routine;
                // 2) The remaining TAUs are set to ZERO by the
                // block factorization routine.
                //
                k = ioffset + jbf;
                //
                // Set INFO on the first occurrence of NaN, NaN takes
                // prcedence over Inf.
                //
                if (iinfo <= n_sub && iinfo > 0) {
                    info = ioffset + iinfo;
                }
                //
                // Return from the routine.
                //
                work[1 - 1] = COMPLEX(lwkopt);
                //
                return;
                //
            }
            //
            j += jbf;
            //
        }
        //
    }
    //
    // Use unblocked code to factor the last or only block.
    // J = JMAX+1 means we factorized the maximum possible number of
    // columns, that is in ELSE clause we need to compute
    // the MAXC2NORM and RELMAXC2NORM to return after we processed
    // the blocks.
    //
    INTEGER kf = 0;
    INTEGER jmaxc2nrm = 0;
    if (j <= jmax) {
        //
        // N_SUB is the number of columns in the submatrix;
        // IOFFSET is the number of rows that should not be factorized.
        //
        n_sub = n - j + 1;
        ioffset = j - 1;
        //
        Claqp2rk(m, n_sub, nrhs, ioffset, jmax - j + 1, abstol, reltol, kp1, maxc2nrm, &a[(j - 1) * lda], lda, kf, maxc2nrmk, relmaxc2nrmk, &jpiv[j - 1], &tau[j - 1], &rwork[j - 1], &rwork[(n + j) - 1], &work[1 - 1], iinfo);
        //
        // ABSTOL or RELTOL criterion is satisfied when the number of
        // the factorized columns KF is smaller then the  number
        // of columns JMAX-J+1 supplied to be factorized by the
        // unblocked routine, we can return from
        // the routine. Perform the following before returning:
        // a) Set the number of factorized columns K,
        // b) MAXC2NRMK and RELMAXC2NRMK are returned by the
        // unblocked factorization routine above.
        //
        k = j - 1 + kf;
        //
        // Set INFO on the first exception occurence.
        //
        // Set INFO on the first exception occurence of Inf or NaN,
        // (NaN takes precedence over Inf).
        //
        if (iinfo > n_sub && info == 0) {
            info = 2 * ioffset + iinfo;
        } else if (iinfo <= n_sub && iinfo > 0) {
            info = ioffset + iinfo;
        }
        //
    } else {
        //
        // Compute the return values for blocked code.
        //
        // Set the number of factorized columns if the unblocked routine
        // was not called.
        //
        k = jmax;
        //
        // If there exits a residual matrix after the blocked code:
        // 1) compute the values of MAXC2NRMK, RELMAXC2NRMK of the
        // residual matrix, otherwise set them to ZERO;
        // 2) Set TAU(K+1:MINMN) to ZERO.
        //
        if (k < minmn) {
            jmaxc2nrm = k + iRamax(n - k, &rwork[(k + 1) - 1], 1);
            maxc2nrmk = rwork[jmaxc2nrm - 1];
            if (k == 0) {
                relmaxc2nrmk = one;
            } else {
                relmaxc2nrmk = maxc2nrmk / maxc2nrm;
            }
            //
            for (j = k + 1; j <= minmn; j = j + 1) {
                tau[j - 1] = czero;
            }
            //
        } else {
            maxc2nrmk = zero;
            relmaxc2nrmk = zero;
            //
        }
        //
        // END IF( J.LE.JMAX ) THEN
        //
    }
    //
    work[1 - 1] = COMPLEX(lwkopt);
    //
    // End of Cgeqp3rk
    //
}

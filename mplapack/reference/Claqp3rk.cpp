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

// Derived from LAPACK routine ZLAQP3RK.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Claqp3rk(INTEGER const m, INTEGER const n, INTEGER const nrhs, INTEGER const ioffset, INTEGER &nb, REAL const abstol, REAL const reltol, INTEGER const kp1, REAL const maxc2nrm, COMPLEX *a, INTEGER const lda, bool &done, INTEGER &kb, REAL &maxc2nrmk, REAL &relmaxc2nrmk, INTEGER *jpiv, COMPLEX *tau, REAL *vn1, REAL *vn2, COMPLEX *auxv, COMPLEX *f, INTEGER const ldf, INTEGER *iwork, INTEGER &info) {
    //
    // Initialize INFO
    //
    info = 0;
    //
    // MINMNFACT in the smallest dimension of the submatrix
    // A(IOFFSET+1:M,1:N) to be factorized.
    //
    INTEGER minmnfact = min(m - ioffset, n);
    INTEGER minmnupdt = min(m - ioffset, n + nrhs);
    nb = min(nb, minmnfact);
    REAL tol3z = sqrt(Rlamch("Epsilon"));
    REAL hugeval = Rlamch("Overflow");
    //
    // Compute factorization in a while loop over NB columns,
    // K is the column index in the block A(1:M,1:N).
    //
    INTEGER k = 0;
    INTEGER lsticc = 0;
    done = false;
    //
    INTEGER i = 0;
    INTEGER kp = 0;
    INTEGER identifier_if = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const REAL zero = 0.0;
    INTEGER j = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER itemp = 0;
    REAL taunan = 0.0;
    COMPLEX aik = 0.0;
    REAL temp = 0.0;
    const REAL one = 1.0;
    REAL temp2 = 0.0;
    while (k < nb && lsticc == 0) {
        k++;
        i = ioffset + k;
        //
        if (i == 1) {
            //
            // We are at the first column of the original whole matrix A_orig,
            // therefore we use the computed KP1 and MAXC2NRM from the
            // main routine.
            //
            kp = kp1;
            //
        } else {
            //
            // Determine the pivot column in K-th step, i.e. the index
            // of the column with the maximum 2-norm in the
            // submatrix A(I:M,K:N).
            //
            kp = (k - 1) + iRamax(n - k + 1, &vn1[k - 1], 1);
            //
            // Determine the maximum column 2-norm and the relative maximum
            // column 2-norm of the submatrix A(I:M,K:N) in step K.
            //
            maxc2nrmk = vn1[kp - 1];
            //
            // ============================================================
            //
            // Check if the submatrix A(I:M,K:N) contains NaN, set
            // INFO parameter to the column number, where the first NaN
            // is found and return from the routine.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.
            //
            if (Risnan(maxc2nrmk)) {
                //
                done = true;
                //
                // Set KB, the number of factorized partial columns
                // that are non-zero in each step in the block,
                // i.e. the rank of the factor R.
                // Set IF, the number of processed rows in the block, which
                // is the same as the number of processed rows in
                // the original whole matrix A_orig.
                //
                kb = k - 1;
                identifier_if = i - 1;
                info = kb + kp;
                //
                // Set RELMAXC2NRMK to NaN.
                //
                relmaxc2nrmk = maxc2nrmk;
                //
                // There is no need to apply the block reflector to the
                // residual of the matrix A stored in A(KB+1:M,KB+1:N),
                // since the submatrix contains NaN and we stop
                // the computation.
                // But, we need to apply the block reflector to the residual
                // right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the
                // residual right hand sides exist.  This occurs
                // when ( NRHS != 0 AND KB <= (M-IOFFSET) ):
                //
                // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
                // A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H.
                //
                if (nrhs > 0 && kb < (m - ioffset)) {
                    Cgemm("No transpose", "Conjugate transpose", m - identifier_if, nrhs, kb, -cone, &a[((identifier_if + 1) - 1)], lda, &f[((n + 1) - 1)], ldf, cone, &a[((identifier_if + 1) - 1) + ((n + 1) - 1) * lda], lda);
                }
                //
                // There is no need to recompute the 2-norm of the
                // difficult columns, since we stop the factorization.
                //
                // Array TAU(KF+1:MINMNFACT) is not set and contains
                // undefined elements.
                //
                // Return from the routine.
                //
                return;
            }
            //
            // Quick return, if the submatrix A(I:M,K:N) is
            // a zero matrix. We need to check it only if the column index
            // (same as row index) is larger than 1, since the condition
            // for the whole original matrix A_orig is checked in the main
            // routine.
            //
            if (maxc2nrmk == zero) {
                //
                done = true;
                //
                // Set KB, the number of factorized partial columns
                // that are non-zero in each step in the block,
                // i.e. the rank of the factor R.
                // Set IF, the number of processed rows in the block, which
                // is the same as the number of processed rows in
                // the original whole matrix A_orig.
                //
                kb = k - 1;
                identifier_if = i - 1;
                relmaxc2nrmk = zero;
                //
                // There is no need to apply the block reflector to the
                // residual of the matrix A stored in A(KB+1:M,KB+1:N),
                // since the submatrix is zero and we stop the computation.
                // But, we need to apply the block reflector to the residual
                // right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the
                // residual right hand sides exist.  This occurs
                // when ( NRHS != 0 AND KB <= (M-IOFFSET) ):
                //
                // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
                // A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H.
                //
                if (nrhs > 0 && kb < (m - ioffset)) {
                    Cgemm("No transpose", "Conjugate transpose", m - identifier_if, nrhs, kb, -cone, &a[((identifier_if + 1) - 1)], lda, &f[((n + 1) - 1)], ldf, cone, &a[((identifier_if + 1) - 1) + ((n + 1) - 1) * lda], lda);
                }
                //
                // There is no need to recompute the 2-norm of the
                // difficult columns, since we stop the factorization.
                //
                // Set TAUs corresponding to the columns that were not
                // factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = CZERO,
                // which is equivalent to seting TAU(K:MINMNFACT) = CZERO.
                //
                for (j = k; j <= minmnfact; j = j + 1) {
                    tau[j - 1] = czero;
                }
                //
                // Return from the routine.
                //
                return;
                //
            }
            //
            // ============================================================
            //
            // Check if the submatrix A(I:M,K:N) contains Inf,
            // set INFO parameter to the column number, where
            // the first Inf is found plus N, and continue
            // the computation.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.
            //
            if (info == 0 && maxc2nrmk > hugeval) {
                info = n + k - 1 + kp;
            }
            //
            // ============================================================
            //
            // Test for the second and third tolerance stopping criteria.
            // NOTE: There is no need to test for ABSTOL.GE.ZERO, since
            // MAXC2NRMK is non-negative. Similarly, there is no need
            // to test for RELTOL.GE.ZERO, since RELMAXC2NRMK is
            // non-negative.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.
            //
            relmaxc2nrmk = maxc2nrmk / maxc2nrm;
            //
            if (maxc2nrmk <= abstol || relmaxc2nrmk <= reltol) {
                //
                done = true;
                //
                // Set KB, the number of factorized partial columns
                // that are non-zero in each step in the block,
                // i.e. the rank of the factor R.
                // Set IF, the number of processed rows in the block, which
                // is the same as the number of processed rows in
                // the original whole matrix A_orig;
                //
                kb = k - 1;
                identifier_if = i - 1;
                //
                // Apply the block reflector to the residual of the
                // matrix A and the residual of the right hand sides B, if
                // the residual matrix and and/or the residual of the right
                // hand sides exist,  i.e. if the submatrix
                // A(I+1:M,KB+1:N+NRHS) exists.  This occurs when
                // KB < MINMNUPDT = min( M-IOFFSET, N+NRHS ):
                //
                // A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) -
                // A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**H.
                //
                if (kb < minmnupdt) {
                    Cgemm("No transpose", "Conjugate transpose", m - identifier_if, n + nrhs - kb, kb, -cone, &a[((identifier_if + 1) - 1)], lda, &f[((kb + 1) - 1)], ldf, cone, &a[((identifier_if + 1) - 1) + ((kb + 1) - 1) * lda], lda);
                }
                //
                // There is no need to recompute the 2-norm of the
                // difficult columns, since we stop the factorization.
                //
                // Set TAUs corresponding to the columns that were not
                // factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = CZERO,
                // which is equivalent to seting TAU(K:MINMNFACT) = CZERO.
                //
                for (j = k; j <= minmnfact; j = j + 1) {
                    tau[j - 1] = czero;
                }
                //
                // Return from the routine.
                //
                return;
                //
            }
            //
            // ============================================================
            //
            // End ELSE of IF(I.EQ.1)
            //
        }
        //
        // ===============================================================
        //
        // If the pivot column is not the first column of the
        // subblock A(1:M,K:N):
        // 1) swap the K-th column and the KP-th pivot column
        // in A(1:M,1:N);
        // 2) swap the K-th row and the KP-th row in F(1:N,1:K-1)
        // 3) copy the K-th element into the KP-th element of the partial
        // and exact 2-norm vectors VN1 and VN2. (Swap is not needed
        // for VN1 and VN2 since we use the element with the index
        // larger than K in the next loop step.)
        // 4) Save the pivot interchange with the indices relative to the
        // the original matrix A_orig, not the block A(1:M,1:N).
        //
        if (kp != k) {
            Cswap(m, &a[(kp - 1) * lda], 1, &a[(k - 1) * lda], 1);
            Cswap(k - 1, &f[(kp - 1)], ldf, &f[(k - 1)], ldf);
            vn1[kp - 1] = vn1[k - 1];
            vn2[kp - 1] = vn2[k - 1];
            itemp = jpiv[kp - 1];
            jpiv[kp - 1] = jpiv[k - 1];
            jpiv[k - 1] = itemp;
        }
        //
        // Apply previous Householder reflectors to column K:
        // A(I:M,K) := A(I:M,K) - A(I:M,1:K-1)*F(K,1:K-1)**H.
        //
        if (k > 1) {
            for (j = 1; j <= k - 1; j = j + 1) {
                f[(k - 1) + (j - 1) * ldf] = conj(f[(k - 1) + (j - 1) * ldf]);
            }
            Cgemv("No transpose", m - i + 1, k - 1, -cone, &a[(i - 1)], lda, &f[(k - 1)], ldf, cone, &a[(i - 1) + (k - 1) * lda], 1);
            for (j = 1; j <= k - 1; j = j + 1) {
                f[(k - 1) + (j - 1) * ldf] = conj(f[(k - 1) + (j - 1) * ldf]);
            }
        }
        //
        // Generate elementary reflector H(k) using the column A(I:M,K).
        //
        if (i < m) {
            Clarfg(m - i + 1, a[(i - 1) + (k - 1) * lda], &a[((i + 1) - 1) + (k - 1) * lda], 1, tau[k - 1]);
        } else {
            tau[k - 1] = czero;
        }
        //
        // Check if TAU(K) contains NaN, set INFO parameter
        // to the column number where NaN is found and return from
        // the routine.
        // NOTE: There is no need to check TAU(K) for Inf,
        // since Clarfg cannot produce TAU(KK) or Householder vector
        // below the diagonal containing Inf. Only BETA on the diagonal,
        // returned by Clarfg can contain Inf, which requires
        // TAU(K) to contain NaN. Therefore, this case of generating Inf
        // by Clarfg is covered by checking TAU(K) for NaN.
        //
        if (Risnan(tau[k - 1].real())) {
            taunan = tau[k - 1].real();
        } else if (Risnan(tau[k - 1].imag())) {
            taunan = tau[k - 1].imag();
        } else {
            taunan = zero;
        }
        //
        if (Risnan(taunan)) {
            //
            done = true;
            //
            // Set KB, the number of factorized partial columns
            // that are non-zero in each step in the block,
            // i.e. the rank of the factor R.
            // Set IF, the number of processed rows in the block, which
            // is the same as the number of processed rows in
            // the original whole matrix A_orig.
            //
            kb = k - 1;
            identifier_if = i - 1;
            info = k;
            //
            // Set MAXC2NRMK and  RELMAXC2NRMK to NaN.
            //
            maxc2nrmk = taunan;
            relmaxc2nrmk = taunan;
            //
            // There is no need to apply the block reflector to the
            // residual of the matrix A stored in A(KB+1:M,KB+1:N),
            // since the submatrix contains NaN and we stop
            // the computation.
            // But, we need to apply the block reflector to the residual
            // right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the
            // residual right hand sides exist.  This occurs
            // when ( NRHS != 0 AND KB <= (M-IOFFSET) ):
            //
            // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
            // A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H.
            //
            if (nrhs > 0 && kb < (m - ioffset)) {
                Cgemm("No transpose", "Conjugate transpose", m - identifier_if, nrhs, kb, -cone, &a[((identifier_if + 1) - 1)], lda, &f[((n + 1) - 1)], ldf, cone, &a[((identifier_if + 1) - 1) + ((n + 1) - 1) * lda], lda);
            }
            //
            // There is no need to recompute the 2-norm of the
            // difficult columns, since we stop the factorization.
            //
            // Array TAU(KF+1:MINMNFACT) is not set and contains
            // undefined elements.
            //
            // Return from the routine.
            //
            return;
        }
        //
        // ===============================================================
        //
        aik = a[(i - 1) + (k - 1) * lda];
        a[(i - 1) + (k - 1) * lda] = cone;
        //
        // ===============================================================
        //
        // Compute the current K-th column of F:
        // 1) F(K+1:N,K) := tau(K) * A(I:M,K+1:N)**H * A(I:M,K).
        //
        if (k < n + nrhs) {
            Cgemv("Conjugate transpose", m - i + 1, n + nrhs - k, tau[k - 1], &a[(i - 1) + ((k + 1) - 1) * lda], lda, &a[(i - 1) + (k - 1) * lda], 1, czero, &f[((k + 1) - 1) + (k - 1) * ldf], 1);
        }
        //
        // 2) Zero out elements above and on the diagonal of the
        // column K in matrix F, i.e elements F(1:K,K).
        //
        for (j = 1; j <= k; j = j + 1) {
            f[(j - 1) + (k - 1) * ldf] = czero;
        }
        //
        // 3) Incremental updating of the K-th column of F:
        // F(1:N,K) := F(1:N,K) - tau(K) * F(1:N,1:K-1) * A(I:M,1:K-1)**H
        // * A(I:M,K).
        //
        if (k > 1) {
            Cgemv("Conjugate Transpose", m - i + 1, k - 1, -tau[k - 1], &a[(i - 1)], lda, &a[(i - 1) + (k - 1) * lda], 1, czero, &auxv[1 - 1], 1);
            //
            Cgemv("No transpose", n + nrhs, k - 1, cone, &f[0], ldf, &auxv[1 - 1], 1, cone, &f[(k - 1) * ldf], 1);
        }
        //
        // ===============================================================
        //
        // Update the current I-th row of A:
        // A(I,K+1:N+NRHS) := A(I,K+1:N+NRHS)
        // - A(I,1:K)*F(K+1:N+NRHS,1:K)**H.
        //
        if (k < n + nrhs) {
            Cgemm("No transpose", "Conjugate transpose", 1, n + nrhs - k, k, -cone, &a[(i - 1)], lda, &f[((k + 1) - 1)], ldf, cone, &a[(i - 1) + ((k + 1) - 1) * lda], lda);
        }
        //
        a[(i - 1) + (k - 1) * lda] = aik;
        //
        // Update the partial column 2-norms for the residual matrix,
        // only if the residual matrix A(I+1:M,K+1:N) exists, i.e.
        // when K < MINMNFACT = min( M-IOFFSET, N ).
        //
        if (k < minmnfact) {
            //
            for (j = k + 1; j <= n; j = j + 1) {
                if (vn1[j - 1] != zero) {
                    //
                    // NOTE: The following lines follow from the analysis in
                    // Lapack Working Note 176.
                    //
                    temp = abs(a[(i - 1) + (j - 1) * lda]) / vn1[j - 1];
                    temp = max(zero, (one + temp) * (one - temp));
                    temp2 = temp * pow2((vn1[j - 1] / vn2[j - 1]));
                    if (temp2 <= tol3z) {
                        //
                        // At J-index, we have a difficult column for the
                        // update of the 2-norm. Save the index of the previous
                        // difficult column in IWORK(J-1).
                        // NOTE: ILSTCC > 1, threfore we can use IWORK only
                        // with N-1 elements, where the elements are
                        // shifted by 1 to the left.
                        //
                        iwork[(j - 1) - 1] = lsticc;
                        //
                        // Set the index of the last difficult column LSTICC.
                        //
                        lsticc = j;
                        //
                    } else {
                        vn1[j - 1] = vn1[j - 1] * sqrt(temp);
                    }
                }
            }
            //
        }
        //
        // End of while loop.
        //
    }
    //
    // Now, afler the loop:
    // Set KB, the number of factorized columns in the block;
    // Set IF, the number of processed rows in the block, which
    // is the same as the number of processed rows in
    // the original whole matrix A_orig, IF = IOFFSET + KB.
    //
    kb = k;
    identifier_if = i;
    //
    // Apply the block reflector to the residual of the matrix A
    // and the residual of the right hand sides B, if the residual
    // matrix and and/or the residual of the right hand sides
    // exist,  i.e. if the submatrix A(I+1:M,KB+1:N+NRHS) exists.
    // This occurs when KB < MINMNUPDT = min( M-IOFFSET, N+NRHS ):
    //
    // A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) -
    // A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**H.
    //
    if (kb < minmnupdt) {
        Cgemm("No transpose", "Conjugate transpose", m - identifier_if, n + nrhs - kb, kb, -cone, &a[((identifier_if + 1) - 1)], lda, &f[((kb + 1) - 1)], ldf, cone, &a[((identifier_if + 1) - 1) + ((kb + 1) - 1) * lda], lda);
    }
    //
    // Recompute the 2-norm of the difficult columns.
    // Loop over the index of the difficult columns from the largest
    // to the smallest index.
    //
    while (lsticc > 0) {
        //
        // LSTICC is the index of the last difficult column is greater
        // than 1.
        // ITEMP is the index of the previous difficult column.
        //
        itemp = iwork[(lsticc - 1) - 1];
        //
        // Compute the 2-norm explicilty for the last difficult column and
        // save it in the partial and exact 2-norm vectors VN1 and VN2.
        //
        // NOTE: The computation of VN1( LSTICC ) relies on the fact that
        // RCnrm2 does not fail on vectors with norm below the value of
        // SQRT(Rlamch('S'))
        //
        vn1[lsticc - 1] = RCnrm2(m - identifier_if, &a[((identifier_if + 1) - 1) + (lsticc - 1) * lda], 1);
        vn2[lsticc - 1] = vn1[lsticc - 1];
        //
        // Downdate the index of the last difficult column to
        // the index of the previous difficult column.
        //
        lsticc = itemp;
        //
    }
    //
    // End of Claqp3rk
    //
}

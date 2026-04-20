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

// Derived from LAPACK routine ZLAQP2RK.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Claqp2rk(INTEGER const m, INTEGER const n, INTEGER const nrhs, INTEGER const ioffset, INTEGER const kmax, REAL const abstol, REAL const reltol, INTEGER const kp1, REAL const maxc2nrm, COMPLEX *a, INTEGER const lda, INTEGER &k, REAL &maxc2nrmk, REAL &relmaxc2nrmk, INTEGER *jpiv, COMPLEX *tau, REAL *vn1, REAL *vn2, COMPLEX *work, INTEGER &info) {
    //
    // Initialize INFO
    //
    info = 0;
    //
    // MINMNFACT in the smallest dimension of the submatrix
    // A(IOFFSET+1:M,1:N) to be factorized.
    //
    // MINMNUPDT is the smallest dimension
    // of the subarray A(IOFFSET+1:M,1:N+NRHS) to be udated, which
    // contains the submatrices A(IOFFSET+1:M,1:N) and
    // B(IOFFSET+1:M,1:NRHS) as column blocks.
    //
    INTEGER minmnfact = min(m - ioffset, n);
    INTEGER minmnupdt = min(m - ioffset, n + nrhs);
    INTEGER kbound = min(kmax, minmnfact);
    REAL tol3z = sqrt(Rlamch("Epsilon"));
    REAL hugeval = Rlamch("Overflow");
    //
    // Compute the factorization, KK is the lomn loop index.
    //
    INTEGER kk = 0;
    INTEGER i = 0;
    INTEGER kp = 0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER itemp = 0;
    REAL taunan = 0.0;
    const REAL one = 1.0;
    REAL temp = 0.0;
    REAL temp2 = 0.0;
    for (kk = 1; kk <= kbound; kk = kk + 1) {
        //
        i = ioffset + kk;
        //
        if (i == 1) {
            //
            // ============================================================
            //
            // We are at the first column of the original whole matrix A,
            // therefore we use the computed KP1 and MAXC2NRM from the
            // main routine.
            //
            kp = kp1;
            //
            // ============================================================
            //
        } else {
            //
            // ============================================================
            //
            // Determine the pivot column in KK-th step, i.e. the index
            // of the column with the maximum 2-norm in the
            // submatrix A(I:M,K:N).
            //
            kp = (kk - 1) + iRamax(n - kk + 1, &vn1[kk - 1], 1);
            //
            // Determine the maximum column 2-norm and the relative maximum
            // column 2-norm of the submatrix A(I:M,KK:N) in step KK.
            // RELMAXC2NRMK  will be computed later, after somecondition
            // checks on MAXC2NRMK.
            //
            maxc2nrmk = vn1[kp - 1];
            //
            // ============================================================
            //
            // Check if the submatrix A(I:M,KK:N) contains NaN, and set
            // INFO parameter to the column number, where the first NaN
            // is found and return from the routine.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.
            //
            if (Risnan(maxc2nrmk)) {
                //
                // Set K, the number of factorized columns.
                // that are not zero.
                //
                k = kk - 1;
                info = k + kp;
                //
                // Set RELMAXC2NRMK to NaN.
                //
                relmaxc2nrmk = maxc2nrmk;
                //
                // Array TAU(K+1:MINMNFACT) is not set and contains
                // undefined elements.
                //
                return;
            }
            //
            // ============================================================
            //
            // Quick return, if the submatrix A(I:M,KK:N) is
            // a zero matrix.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.
            //
            if (maxc2nrmk == zero) {
                //
                // Set K, the number of factorized columns.
                // that are not zero.
                //
                k = kk - 1;
                relmaxc2nrmk = zero;
                //
                // Set TAUs corresponding to the columns that were not
                // factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to CZERO.
                //
                for (j = kk; j <= minmnfact; j = j + 1) {
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
            // Check if the submatrix A(I:M,KK:N) contains Inf,
            // set INFO parameter to the column number, where
            // the first Inf is found plus N, and continue
            // the computation.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.
            //
            if (info == 0 && maxc2nrmk > hugeval) {
                info = n + kk - 1 + kp;
            }
            //
            // ============================================================
            //
            // Test for the second and third stopping criteria.
            // NOTE: There is no need to test for ABSTOL >= ZERO, since
            // MAXC2NRMK is non-negative. Similarly, there is no need
            // to test for RELTOL >= ZERO, since RELMAXC2NRMK is
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
                // Set K, the number of factorized columns.
                //
                k = kk - 1;
                //
                // Set TAUs corresponding to the columns that were not
                // factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to CZERO.
                //
                for (j = kk; j <= minmnfact; j = j + 1) {
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
        // subblock A(1:M,KK:N):
        // 1) swap the KK-th column and the KP-th pivot column
        // in A(1:M,1:N);
        // 2) copy the KK-th element into the KP-th element of the partial
        // and exact 2-norm vectors VN1 and VN2. ( Swap is not needed
        // for VN1 and VN2 since we use the element with the index
        // larger than KK in the next loop step.)
        // 3) Save the pivot interchange with the indices relative to the
        // the original matrix A, not the block A(1:M,1:N).
        //
        if (kp != kk) {
            Cswap(m, &a[(kp - 1) * lda], 1, &a[(kk - 1) * lda], 1);
            vn1[kp - 1] = vn1[kk - 1];
            vn2[kp - 1] = vn2[kk - 1];
            itemp = jpiv[kp - 1];
            jpiv[kp - 1] = jpiv[kk - 1];
            jpiv[kk - 1] = itemp;
        }
        //
        // Generate elementary reflector H(KK) using the column A(I:M,KK),
        // if the column has more than one element, otherwise
        // the elementary reflector would be an identity matrix,
        // and TAU(KK) = CZERO.
        //
        if (i < m) {
            Clarfg(m - i + 1, a[(i - 1) + (kk - 1) * lda], &a[((i + 1) - 1) + (kk - 1) * lda], 1, tau[kk - 1]);
        } else {
            tau[kk - 1] = czero;
        }
        //
        // Check if TAU(KK) contains NaN, set INFO parameter
        // to the column number where NaN is found and return from
        // the routine.
        // NOTE: There is no need to check TAU(KK) for Inf,
        // since Clarfg cannot produce TAU(KK) or Householder vector
        // below the diagonal containing Inf. Only BETA on the diagonal,
        // returned by Clarfg can contain Inf, which requires
        // TAU(KK) to contain NaN. Therefore, this case of generating Inf
        // by Clarfg is covered by checking TAU(KK) for NaN.
        //
        if (Risnan(tau[kk - 1].real())) {
            taunan = tau[kk - 1].real();
        } else if (Risnan(tau[kk - 1].imag())) {
            taunan = tau[kk - 1].imag();
        } else {
            taunan = zero;
        }
        //
        if (Risnan(taunan)) {
            k = kk - 1;
            info = kk;
            //
            // Set MAXC2NRMK and  RELMAXC2NRMK to NaN.
            //
            maxc2nrmk = taunan;
            relmaxc2nrmk = taunan;
            //
            // Array TAU(KK:MINMNFACT) is not set and contains
            // undefined elements, except the first element TAU(KK) = NaN.
            //
            return;
        }
        //
        // Apply H(KK)**H to A(I:M,KK+1:N+NRHS) from the left.
        // ( If M >= N, then at KK = N there is no residual matrix,
        // i.e. no columns of A to update, only columns of B.
        // If M < N, then at KK = M-IOFFSET, I = M and we have a
        // one-row residual matrix in A and the elementary
        // reflector is a unit matrix, TAU(KK) = CZERO, i.e. no update
        // is needed for the residual matrix in A and the
        // right-hand-side-matrix in B.
        // Therefore, we update only if
        // KK < MINMNUPDT = min(M-IOFFSET, N+NRHS)
        // condition is satisfied, not only KK < N+NRHS )
        //
        if (kk < minmnupdt) {
            Clarf1f("Left", m - i + 1, n + nrhs - kk, &a[(i - 1) + (kk - 1) * lda], 1, conj(tau[kk - 1]), &a[(i - 1) + ((kk + 1) - 1) * lda], lda, &work[1 - 1]);
        }
        //
        if (kk < minmnfact) {
            //
            // Update the partial column 2-norms for the residual matrix,
            // only if the residual matrix A(I+1:M,KK+1:N) exists, i.e.
            // when KK < min(M-IOFFSET, N).
            //
            for (j = kk + 1; j <= n; j = j + 1) {
                if (vn1[j - 1] != zero) {
                    //
                    // NOTE: The following lines follow from the analysis in
                    // Lapack Working Note 176.
                    //
                    temp = one - pow2((abs(a[(i - 1) + (j - 1) * lda]) / vn1[j - 1]));
                    temp = max(temp, zero);
                    temp2 = temp * pow2((vn1[j - 1] / vn2[j - 1]));
                    if (temp2 <= tol3z) {
                        //
                        // Compute the column 2-norm for the partial
                        // column A(I+1:M,J) by explicitly computing it,
                        // and store it in both partial 2-norm vector VN1
                        // and exact column 2-norm vector VN2.
                        //
                        vn1[j - 1] = RCnrm2(m - i, &a[((i + 1) - 1) + (j - 1) * lda], 1);
                        vn2[j - 1] = vn1[j - 1];
                        //
                    } else {
                        //
                        // Update the column 2-norm for the partial
                        // column A(I+1:M,J) by removing one
                        // element A(I,J) and store it in partial
                        // 2-norm vector VN1.
                        //
                        vn1[j - 1] = vn1[j - 1] * sqrt(temp);
                        //
                    }
                }
            }
            //
        }
        //
        // End factorization loop
        //
    }
    //
    // If we reached this point, all colunms have been factorized,
    // i.e. no condition was triggered to exit the routine.
    // Set the number of factorized columns.
    //
    k = kbound;
    //
    // We reached the end of the loop, i.e. all KMAX columns were
    // factorized, we need to set MAXC2NRMK and RELMAXC2NRMK before
    // we return.
    //
    INTEGER jmaxc2nrm = 0;
    if (k < minmnfact) {
        //
        jmaxc2nrm = k + iRamax(n - k, &vn1[(k + 1) - 1], 1);
        maxc2nrmk = vn1[jmaxc2nrm - 1];
        //
        if (k == 0) {
            relmaxc2nrmk = one;
        } else {
            relmaxc2nrmk = maxc2nrmk / maxc2nrm;
        }
        //
    } else {
        maxc2nrmk = zero;
        relmaxc2nrmk = zero;
    }
    //
    // We reached the end of the loop, i.e. all KMAX columns were
    // factorized, set TAUs corresponding to the columns that were
    // not factorized to ZERO, i.e. TAU(K+1:MINMNFACT) set to CZERO.
    //
    for (j = k + 1; j <= minmnfact; j = j + 1) {
        tau[j - 1] = czero;
    }
    //
    // End of Claqp2rk
    //
}

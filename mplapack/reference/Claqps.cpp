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

void Claqps(INTEGER const m, INTEGER const n, INTEGER const offset, INTEGER const nb, INTEGER &kb, COMPLEX *a, INTEGER const lda, INTEGER *jpvt, COMPLEX *tau, REAL *vn1, REAL *vn2, COMPLEX *auxv, COMPLEX *f, INTEGER const ldf) {
    INTEGER lastrk = 0;
    INTEGER lsticc = 0;
    INTEGER k = 0;
    REAL tol3z = 0.0;
    INTEGER rk = 0;
    INTEGER pvt = 0;
    INTEGER itemp = 0;
    INTEGER j = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    COMPLEX akk = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const REAL zero = 0.0;
    REAL temp = 0.0;
    const REAL one = 1.0;
    REAL temp2 = 0.0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    lastrk = min(m, n + offset);
    lsticc = 0;
    k = 0;
    tol3z = sqrt(Rlamch("Epsilon"));
//
//     Beginning of while loop.
//
statement_10:
    if ((k < nb) && (lsticc == 0)) {
        k++;
        rk = offset + k;
        //
        //        Determine ith pivot column and swap if necessary
        //
        pvt = (k - 1) + iRamax(n - k + 1, &vn1[k - 1], 1);
        if (pvt != k) {
            Cswap(m, &a[(pvt - 1) * lda], 1, &a[(k - 1) * lda], 1);
            Cswap(k - 1, &f[(pvt - 1)], ldf, &f[(k - 1)], ldf);
            itemp = jpvt[pvt - 1];
            jpvt[pvt - 1] = jpvt[k - 1];
            jpvt[k - 1] = itemp;
            vn1[pvt - 1] = vn1[k - 1];
            vn2[pvt - 1] = vn2[k - 1];
        }
        //
        //        Apply previous Householder reflectors to column K:
        //        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**H.
        //
        if (k > 1) {
            for (j = 1; j <= k - 1; j = j + 1) {
                f[(k - 1) + (j - 1) * ldf] = conj(f[(k - 1) + (j - 1) * ldf]);
            }
            Cgemv("No transpose", m - rk + 1, k - 1, -cone, &a[(rk - 1)], lda, &f[(k - 1)], ldf, cone, &a[(rk - 1) + (k - 1) * lda], 1);
            for (j = 1; j <= k - 1; j = j + 1) {
                f[(k - 1) + (j - 1) * ldf] = conj(f[(k - 1) + (j - 1) * ldf]);
            }
        }
        //
        //        Generate elementary reflector H(k).
        //
        if (rk < m) {
            Clarfg(m - rk + 1, a[(rk - 1) + (k - 1) * lda], &a[((rk + 1) - 1) + (k - 1) * lda], 1, tau[k - 1]);
        } else {
            Clarfg(1, a[(rk - 1) + (k - 1) * lda], &a[(rk - 1) + (k - 1) * lda], 1, tau[k - 1]);
        }
        //
        akk = a[(rk - 1) + (k - 1) * lda];
        a[(rk - 1) + (k - 1) * lda] = cone;
        //
        //        Compute Kth column of F:
        //
        //        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).
        //
        if (k < n) {
            Cgemv("Conjugate transpose", m - rk + 1, n - k, tau[k - 1], &a[(rk - 1) + ((k + 1) - 1) * lda], lda, &a[(rk - 1) + (k - 1) * lda], 1, czero, &f[((k + 1) - 1) + (k - 1) * ldf], 1);
        }
        //
        //        Padding F(1:K,K) with zeros.
        //
        for (j = 1; j <= k; j = j + 1) {
            f[(j - 1) + (k - 1) * ldf] = czero;
        }
        //
        //        Incremental updating of F:
        //        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**H
        //                    *A(RK:M,K).
        //
        if (k > 1) {
            Cgemv("Conjugate transpose", m - rk + 1, k - 1, -tau[k - 1], &a[(rk - 1)], lda, &a[(rk - 1) + (k - 1) * lda], 1, czero, &auxv[1 - 1], 1);
            //
            Cgemv("No transpose", n, k - 1, cone, &f[(1 - 1)], ldf, &auxv[1 - 1], 1, cone, &f[(k - 1) * ldf], 1);
        }
        //
        //        Update the current row of A:
        //        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H.
        //
        if (k < n) {
            Cgemm("No transpose", "Conjugate transpose", 1, n - k, k, -cone, &a[(rk - 1)], lda, &f[((k + 1) - 1)], ldf, cone, &a[(rk - 1) + ((k + 1) - 1) * lda], lda);
        }
        //
        //        Update partial column norms.
        //
        if (rk < lastrk) {
            for (j = k + 1; j <= n; j = j + 1) {
                if (vn1[j - 1] != zero) {
                    //
                    //                 NOTE: The following 4 lines follow from the analysis in
                    //                 Lapack Working Note 176.
                    //
                    temp = abs(a[(rk - 1) + (j - 1) * lda]) / vn1[j - 1];
                    temp = max(zero, (one + temp) * (one - temp));
                    temp2 = temp * pow2((vn1[j - 1] / vn2[j - 1]));
                    if (temp2 <= tol3z) {
                        vn2[j - 1] = castREAL(lsticc);
                        lsticc = j;
                    } else {
                        vn1[j - 1] = vn1[j - 1] * sqrt(temp);
                    }
                }
            }
        }
        //
        a[(rk - 1) + (k - 1) * lda] = akk;
        //
        //        End of while loop.
        //
        goto statement_10;
    }
    kb = k;
    rk = offset + kb;
    //
    //     Apply the block reflector to the rest of the matrix:
    //     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
    //                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**H.
    //
    if (kb < min(n, m - offset)) {
        Cgemm("No transpose", "Conjugate transpose", m - rk, n - kb, kb, -cone, &a[((rk + 1) - 1)], lda, &f[((kb + 1) - 1)], ldf, cone, &a[((rk + 1) - 1) + ((kb + 1) - 1) * lda], lda);
    }
//
//     Recomputation of difficult columns.
//
statement_60:
    if (lsticc > 0) {
        itemp = nint(vn2[lsticc - 1]);
        vn1[lsticc - 1] = RCnrm2(m - rk, &a[((rk + 1) - 1) + (lsticc - 1) * lda], 1);
        //
        //        NOTE: The computation of VN1( LSTICC ) relies on the fact that
        //        SNRM2 does not fail on vectors with norm below the value of
        //        SQRT(DLAMCH('S'))
        //
        vn2[lsticc - 1] = vn1[lsticc - 1];
        lsticc = itemp;
        goto statement_60;
    }
    //
    //     End of Claqps
    //
}

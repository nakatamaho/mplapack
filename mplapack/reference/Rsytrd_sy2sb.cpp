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

void Rsytrd_sy2sb(const char *uplo, INTEGER const n, INTEGER const kd, REAL *a, INTEGER const lda, REAL *ab, INTEGER const ldab, REAL *tau, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Determine the minimal workspace size required
    //     and test the input parameters
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    bool lquery = (lwork == -1);
    INTEGER lwmin = iMlaenv2stage(4, "Rsytrd_sy2sb", "", n, kd, -1, -1);
    //
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kd < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldab < max((INTEGER)1, kd + 1)) {
        info = -7;
    } else if (lwork < lwmin && !lquery) {
        info = -10;
    }
    //
    if (info != 0) {
        Mxerbla("Rsytrd_sy2sb", -info);
        return;
    } else if (lquery) {
        work[1 - 1] = lwmin;
        return;
    }
    //
    //     Quick return if possible
    //     Copy the upper/lower portion of A into AB
    //
    INTEGER i = 0;
    INTEGER lk = 0;
    if (n <= kd + 1) {
        if (upper) {
            for (i = 1; i <= n; i = i + 1) {
                lk = min(kd + 1, i);
                Rcopy(lk, &a[((i - lk + 1) - 1) + (i - 1) * lda], 1, &ab[((kd + 1 - lk + 1) - 1) + (i - 1) * ldab], 1);
            }
        } else {
            for (i = 1; i <= n; i = i + 1) {
                lk = min(kd + 1, n - i + 1);
                Rcopy(lk, &a[(i - 1) + (i - 1) * lda], 1, &ab[(i - 1) * ldab], 1);
            }
        }
        work[1 - 1] = 1;
        return;
    }
    //
    //     Determine the pointer position for the workspace
    //
    INTEGER ldt = kd;
    INTEGER lds1 = kd;
    INTEGER lt = ldt * kd;
    INTEGER lw = n * kd;
    INTEGER ls1 = lds1 * kd;
    INTEGER ls2 = lwmin - lt - lw - ls1;
    //      LS2 = N*MAX(KD,FACTOPTNB)
    INTEGER tpos = 1;
    INTEGER wpos = tpos + lt;
    INTEGER s1pos = wpos + lw;
    INTEGER s2pos = s1pos + ls1;
    INTEGER ldw = 0;
    INTEGER lds2 = 0;
    if (upper) {
        ldw = kd;
        lds2 = kd;
    } else {
        ldw = n;
        lds2 = n;
    }
    //
    //     Set the workspace of the triangular matrix T to zero once such a
    //     way every time T is generated the upper/lower portion will be always zero
    //
    const REAL zero = 0.0;
    Rlaset("A", ldt, kd, zero, zero, &work[tpos - 1], ldt);
    //
    INTEGER pn = 0;
    INTEGER pk = 0;
    INTEGER iinfo = 0;
    INTEGER j = 0;
    const REAL one = 1.0;
    const REAL half = 0.5e+0;
    const REAL rone = 1.0;
    if (upper) {
        for (i = 1; i <= n - kd; i = i + kd) {
            pn = n - i - kd + 1;
            pk = min(n - i - kd + 1, kd);
            //
            //            Compute the LQ factorization of the current block
            //
            Rgelqf(kd, pn, &a[(i - 1) + ((i + kd) - 1) * lda], lda, &tau[i - 1], &work[s2pos - 1], ls2, iinfo);
            //
            //            Copy the upper portion of A into AB
            //
            for (j = i; j <= i + pk - 1; j = j + 1) {
                lk = min(kd, n - j) + 1;
                Rcopy(lk, &a[(j - 1) + (j - 1) * lda], lda, &ab[((kd + 1) - 1) + (j - 1) * ldab], ldab - 1);
            }
            //
            Rlaset("Lower", pk, pk, zero, one, &a[(i - 1) + ((i + kd) - 1) * lda], lda);
            //
            //            Form the matrix T
            //
            Rlarft("Forward", "Rowwise", pn, pk, &a[(i - 1) + ((i + kd) - 1) * lda], lda, &tau[i - 1], &work[tpos - 1], ldt);
            //
            //            Compute W:
            //
            Rgemm("Conjugate", "No transpose", pk, pn, pk, one, &work[tpos - 1], ldt, &a[(i - 1) + ((i + kd) - 1) * lda], lda, zero, &work[s2pos - 1], lds2);
            //
            Rsymm("Right", uplo, pk, pn, one, &a[((i + kd) - 1) + ((i + kd) - 1) * lda], lda, &work[s2pos - 1], lds2, zero, &work[wpos - 1], ldw);
            //
            Rgemm("No transpose", "Conjugate", pk, pk, pn, one, &work[wpos - 1], ldw, &work[s2pos - 1], lds2, zero, &work[s1pos - 1], lds1);
            //
            Rgemm("No transpose", "No transpose", pk, pn, pk, -half, &work[s1pos - 1], lds1, &a[(i - 1) + ((i + kd) - 1) * lda], lda, one, &work[wpos - 1], ldw);
            //
            //            Update the unreduced submatrix A(i+kd:n,i+kd:n), using
            //            an update of the form:  A := A - V'*W - W'*V
            //
            Rsyr2k(uplo, "Conjugate", pn, pk, -one, &a[(i - 1) + ((i + kd) - 1) * lda], lda, &work[wpos - 1], ldw, rone, &a[((i + kd) - 1) + ((i + kd) - 1) * lda], lda);
        }
        //
        //        Copy the upper band to AB which is the band storage matrix
        //
        for (j = n - kd + 1; j <= n; j = j + 1) {
            lk = min(kd, n - j) + 1;
            Rcopy(lk, &a[(j - 1) + (j - 1) * lda], lda, &ab[((kd + 1) - 1) + (j - 1) * ldab], ldab - 1);
        }
        //
    } else {
        //
        //         Reduce the lower triangle of A to lower band matrix
        //
        for (i = 1; i <= n - kd; i = i + kd) {
            pn = n - i - kd + 1;
            pk = min(n - i - kd + 1, kd);
            //
            //            Compute the QR factorization of the current block
            //
            Rgeqrf(pn, kd, &a[((i + kd) - 1) + (i - 1) * lda], lda, &tau[i - 1], &work[s2pos - 1], ls2, iinfo);
            //
            //            Copy the upper portion of A into AB
            //
            for (j = i; j <= i + pk - 1; j = j + 1) {
                lk = min(kd, n - j) + 1;
                Rcopy(lk, &a[(j - 1) + (j - 1) * lda], 1, &ab[(j - 1) * ldab], 1);
            }
            //
            Rlaset("Upper", pk, pk, zero, one, &a[((i + kd) - 1) + (i - 1) * lda], lda);
            //
            //            Form the matrix T
            //
            Rlarft("Forward", "Columnwise", pn, pk, &a[((i + kd) - 1) + (i - 1) * lda], lda, &tau[i - 1], &work[tpos - 1], ldt);
            //
            //            Compute W:
            //
            Rgemm("No transpose", "No transpose", pn, pk, pk, one, &a[((i + kd) - 1) + (i - 1) * lda], lda, &work[tpos - 1], ldt, zero, &work[s2pos - 1], lds2);
            //
            Rsymm("Left", uplo, pn, pk, one, &a[((i + kd) - 1) + ((i + kd) - 1) * lda], lda, &work[s2pos - 1], lds2, zero, &work[wpos - 1], ldw);
            //
            Rgemm("Conjugate", "No transpose", pk, pk, pn, one, &work[s2pos - 1], lds2, &work[wpos - 1], ldw, zero, &work[s1pos - 1], lds1);
            //
            Rgemm("No transpose", "No transpose", pn, pk, pk, -half, &a[((i + kd) - 1) + (i - 1) * lda], lda, &work[s1pos - 1], lds1, one, &work[wpos - 1], ldw);
            //
            //            Update the unreduced submatrix A(i+kd:n,i+kd:n), using
            //            an update of the form:  A := A - V*W' - W*V'
            //
            Rsyr2k(uplo, "No transpose", pn, pk, -one, &a[((i + kd) - 1) + (i - 1) * lda], lda, &work[wpos - 1], ldw, rone, &a[((i + kd) - 1) + ((i + kd) - 1) * lda], lda);
            //            ==================================================================
            //            RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
            //             DO 45 J = I, I+PK-1
            //                LK = MIN( KD, N-J ) + 1
            //                CALL Rcopy( LK, AB( 1, J ), 1, A( J, J ), 1 )
            //   45        CONTINUE
            //            ==================================================================
        }
        //
        //        Copy the lower band to AB which is the band storage matrix
        //
        for (j = n - kd + 1; j <= n; j = j + 1) {
            lk = min(kd, n - j) + 1;
            Rcopy(lk, &a[(j - 1) + (j - 1) * lda], 1, &ab[(j - 1) * ldab], 1);
        }
        //
    }
    //
    work[1 - 1] = lwmin;
    //
    //     End of Rsytrd_sy2sb
    //
}

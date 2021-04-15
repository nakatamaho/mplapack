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

void Rlagv2(REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *alphar, REAL *alphai, REAL *beta, REAL &csl, REAL &snl, REAL &csr, REAL &snr) {
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL safmin = Rlamch("S");
    REAL ulp = Rlamch("P");
    //
    //     Scale A
    //
    REAL anorm = max({abs(a[(1 - 1)]) + abs(a[(2 - 1)]), abs(a[(2 - 1) * lda]) + abs(a[(2 - 1) + (2 - 1) * lda]), safmin});
    const REAL one = 1.0;
    REAL ascale = one / anorm;
    a[(1 - 1)] = ascale * a[(1 - 1)];
    a[(2 - 1) * lda] = ascale * a[(2 - 1) * lda];
    a[(2 - 1)] = ascale * a[(2 - 1)];
    a[(2 - 1) + (2 - 1) * lda] = ascale * a[(2 - 1) + (2 - 1) * lda];
    //
    //     Scale B
    //
    REAL bnorm = max({abs(b[(1 - 1)]), abs(b[(2 - 1) * ldb]) + abs(b[(2 - 1) + (2 - 1) * ldb]), safmin});
    REAL bscale = one / bnorm;
    b[(1 - 1)] = bscale * b[(1 - 1)];
    b[(2 - 1) * ldb] = bscale * b[(2 - 1) * ldb];
    b[(2 - 1) + (2 - 1) * ldb] = bscale * b[(2 - 1) + (2 - 1) * ldb];
    //
    //     Check if A can be deflated
    //
    const REAL zero = 0.0;
    REAL wi = 0.0;
    REAL r = 0.0;
    REAL t = 0.0;
    REAL scale1 = 0.0;
    REAL scale2 = 0.0;
    REAL wr1 = 0.0;
    REAL wr2 = 0.0;
    REAL h1 = 0.0;
    REAL h2 = 0.0;
    REAL h3 = 0.0;
    REAL rr = 0.0;
    REAL qq = 0.0;
    if (abs(a[(2 - 1)]) <= ulp) {
        csl = one;
        snl = zero;
        csr = one;
        snr = zero;
        a[(2 - 1)] = zero;
        b[(2 - 1)] = zero;
        wi = zero;
        //
        //     Check if B is singular
        //
    } else if (abs(b[(1 - 1)]) <= ulp) {
        Rlartg(a[(1 - 1)], a[(2 - 1)], csl, snl, r);
        csr = one;
        snr = zero;
        Rrot(2, &a[(1 - 1)], lda, &a[(2 - 1)], lda, csl, snl);
        Rrot(2, &b[(1 - 1)], ldb, &b[(2 - 1)], ldb, csl, snl);
        a[(2 - 1)] = zero;
        b[(1 - 1)] = zero;
        b[(2 - 1)] = zero;
        wi = zero;
        //
    } else if (abs(b[(2 - 1) + (2 - 1) * ldb]) <= ulp) {
        Rlartg(a[(2 - 1) + (2 - 1) * lda], a[(2 - 1)], csr, snr, t);
        snr = -snr;
        Rrot(2, &a[(1 - 1)], 1, &a[(2 - 1) * lda], 1, csr, snr);
        Rrot(2, &b[(1 - 1)], 1, &b[(2 - 1) * ldb], 1, csr, snr);
        csl = one;
        snl = zero;
        a[(2 - 1)] = zero;
        b[(2 - 1)] = zero;
        b[(2 - 1) + (2 - 1) * ldb] = zero;
        wi = zero;
        //
    } else {
        //
        //        B is nonsingular, first compute the eigenvalues of (A,B)
        //
        Rlag2(a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi);
        //
        if (wi == zero) {
            //
            //           two real eigenvalues, compute s*A-w*B
            //
            h1 = scale1 * a[(1 - 1)] - wr1 * b[(1 - 1)];
            h2 = scale1 * a[(2 - 1) * lda] - wr1 * b[(2 - 1) * ldb];
            h3 = scale1 * a[(2 - 1) + (2 - 1) * lda] - wr1 * b[(2 - 1) + (2 - 1) * ldb];
            //
            rr = Rlapy2(h1, h2);
            qq = Rlapy2(scale1 * a[(2 - 1)], h3);
            //
            if (rr > qq) {
                //
                //              find right rotation matrix to zero 1,1 element of
                //              (sA - wB)
                //
                Rlartg(h2, h1, csr, snr, t);
                //
            } else {
                //
                //              find right rotation matrix to zero 2,1 element of
                //              (sA - wB)
                //
                Rlartg(h3, scale1 * a[(2 - 1)], csr, snr, t);
                //
            }
            //
            snr = -snr;
            Rrot(2, &a[(1 - 1)], 1, &a[(2 - 1) * lda], 1, csr, snr);
            Rrot(2, &b[(1 - 1)], 1, &b[(2 - 1) * ldb], 1, csr, snr);
            //
            //           compute inf norms of A and B
            //
            h1 = max(abs(a[(1 - 1)]) + abs(a[(2 - 1) * lda]), abs(a[(2 - 1)]) + abs(a[(2 - 1) + (2 - 1) * lda]));
            h2 = max(abs(b[(1 - 1)]) + abs(b[(2 - 1) * ldb]), abs(b[(2 - 1)]) + abs(b[(2 - 1) + (2 - 1) * ldb]));
            //
            if ((scale1 * h1) >= abs(wr1) * h2) {
                //
                //              find left rotation matrix Q to zero out B(2,1)
                //
                Rlartg(b[(1 - 1)], b[(2 - 1)], csl, snl, r);
                //
            } else {
                //
                //              find left rotation matrix Q to zero out A(2,1)
                //
                Rlartg(a[(1 - 1)], a[(2 - 1)], csl, snl, r);
                //
            }
            //
            Rrot(2, &a[(1 - 1)], lda, &a[(2 - 1)], lda, csl, snl);
            Rrot(2, &b[(1 - 1)], ldb, &b[(2 - 1)], ldb, csl, snl);
            //
            a[(2 - 1)] = zero;
            b[(2 - 1)] = zero;
            //
        } else {
            //
            //           a pair of complex conjugate eigenvalues
            //           first compute the SVD of the matrix B
            //
            Rlasv2(b[(1 - 1)], b[(2 - 1) * ldb], b[(2 - 1) + (2 - 1) * ldb], r, t, snr, csr, snl, csl);
            //
            //           Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and
            //           Z is right rotation matrix computed from Rlasv2
            //
            Rrot(2, &a[(1 - 1)], lda, &a[(2 - 1)], lda, csl, snl);
            Rrot(2, &b[(1 - 1)], ldb, &b[(2 - 1)], ldb, csl, snl);
            Rrot(2, &a[(1 - 1)], 1, &a[(2 - 1) * lda], 1, csr, snr);
            Rrot(2, &b[(1 - 1)], 1, &b[(2 - 1) * ldb], 1, csr, snr);
            //
            b[(2 - 1)] = zero;
            b[(2 - 1) * ldb] = zero;
            //
        }
        //
    }
    //
    //     Unscaling
    //
    a[(1 - 1)] = anorm * a[(1 - 1)];
    a[(2 - 1)] = anorm * a[(2 - 1)];
    a[(2 - 1) * lda] = anorm * a[(2 - 1) * lda];
    a[(2 - 1) + (2 - 1) * lda] = anorm * a[(2 - 1) + (2 - 1) * lda];
    b[(1 - 1)] = bnorm * b[(1 - 1)];
    b[(2 - 1)] = bnorm * b[(2 - 1)];
    b[(2 - 1) * ldb] = bnorm * b[(2 - 1) * ldb];
    b[(2 - 1) + (2 - 1) * ldb] = bnorm * b[(2 - 1) + (2 - 1) * ldb];
    //
    if (wi == zero) {
        alphar[1 - 1] = a[(1 - 1)];
        alphar[2 - 1] = a[(2 - 1) + (2 - 1) * lda];
        alphai[1 - 1] = zero;
        alphai[2 - 1] = zero;
        beta[1 - 1] = b[(1 - 1)];
        beta[2 - 1] = b[(2 - 1) + (2 - 1) * ldb];
    } else {
        alphar[1 - 1] = anorm * wr1 / scale1 / bnorm;
        alphai[1 - 1] = anorm * wi / scale1 / bnorm;
        alphar[2 - 1] = alphar[1 - 1];
        alphai[2 - 1] = -alphai[1 - 1];
        beta[1 - 1] = one;
        beta[2 - 1] = one;
    }
    //
    //     End of Rlagv2
    //
}

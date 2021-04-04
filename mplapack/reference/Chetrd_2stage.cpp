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

void Chetrd_2stage(const char *vect, const char *uplo, INTEGER const &n, COMPLEX *a, INTEGER const &lda, REAL *d, REAL *e, COMPLEX *tau, COMPLEX *hous2, INTEGER const &lhous2, COMPLEX *work, INTEGER const &lwork, INTEGER &info) {
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
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    bool wantq = Mlsame(vect, "V");
    bool upper = Mlsame(uplo, "U");
    bool lquery = (lwork == -1) || (lhous2 == -1);
    //
    //     Determine the block size, the workspace size and the hous size.
    //
    INTEGER kd = iMlaenv2stage[("Chetrd_2stage" - 1) * ldiMlaenv2stage];
    INTEGER ib = iMlaenv2stage[(2 - 1) + ("Chetrd_2stage" - 1) * ldiMlaenv2stage];
    INTEGER lhmin = iMlaenv2stage[(3 - 1) + ("Chetrd_2stage" - 1) * ldiMlaenv2stage];
    INTEGER lwmin = iMlaenv2stage[(4 - 1) + ("Chetrd_2stage" - 1) * ldiMlaenv2stage];
    //      WRITE(*,*),'Chetrd_2stage N KD UPLO LHMIN LWMIN ',N, KD, UPLO,
    //     $            LHMIN, LWMIN
    //
    if (!Mlsame(vect, "N")) {
        info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (lhous2 < lhmin && !lquery) {
        info = -10;
    } else if (lwork < lwmin && !lquery) {
        info = -12;
    }
    //
    if (info == 0) {
        hous2[1 - 1] = lhmin;
        work[1 - 1] = lwmin;
    }
    //
    if (info != 0) {
        Mxerbla("Chetrd_2stage", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    //     Determine poINTEGERer position
    //
    INTEGER ldab = kd + 1;
    INTEGER lwrk = lwork - ldab * n;
    INTEGER abpos = 1;
    INTEGER wpos = abpos + ldab * n;
    Chetrd_he2hb(uplo, n, kd, a, lda, work[abpos - 1], ldab, tau, work[wpos - 1], lwrk, info);
    if (info != 0) {
        Mxerbla("Chetrd_he2hb", -info);
        return;
    }
    Chetrd_hb2st("Y", vect, uplo, n, kd, work[abpos - 1], ldab, d, e, hous2, lhous2, work[wpos - 1], lwrk, info);
    if (info != 0) {
        Mxerbla("Chetrd_HB2ST", -info);
        return;
    }
    //
    hous2[1 - 1] = lhmin;
    work[1 - 1] = lwmin;
    //
    //     End of Chetrd_2stage
    //
}

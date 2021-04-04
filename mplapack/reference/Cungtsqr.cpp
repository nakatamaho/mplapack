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

void Cungtsqr(INTEGER const &m, INTEGER const &n, INTEGER const &mb, INTEGER const &nb, COMPLEX *a, INTEGER const &lda, COMPLEX *t, INTEGER const &ldt, COMPLEX *work, INTEGER const &lwork, INTEGER &info) {
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
    //     Test the input parameters
    //
    bool lquery = lwork == -1;
    info = 0;
    INTEGER nblocal = 0;
    INTEGER ldc = 0;
    INTEGER lc = 0;
    INTEGER lw = 0;
    INTEGER lworkopt = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0 || m < n) {
        info = -2;
    } else if (mb <= n) {
        info = -3;
    } else if (nb < 1) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -6;
    } else if (ldt < max((INTEGER)1, min(nb, n))) {
        info = -8;
    } else {
        //
        //        This workspace is used to store array C(LDC, N) and WORK(LWORK)
        //        in the call to Clamtsqr. See the documentation for Clamtsqr.
        //
        if (lwork < 2 && (!lquery)) {
            info = -10;
        } else {
            //
            //           Set block size for column blocks
            //
            nblocal = min(nb, n);
            //
            //           LWORK = -1, then set the size for the array C(LDC,N)
            //           in Clamtsqr call and set the optimal size of the work array
            //           WORK(LWORK) in Clamtsqr call.
            //
            ldc = m;
            lc = ldc * n;
            lw = n * nblocal;
            //
            lworkopt = lc + lw;
            //
            if ((lwork < max((INTEGER)1, lworkopt)) && (!lquery)) {
                info = -10;
            }
        }
        //
    }
    //
    //     Handle error in the input parameters and return workspace query.
    //
    if (info != 0) {
        Mxerbla("Cungtsqr", -info);
        return;
    } else if (lquery) {
        work[1 - 1] = COMPLEX(lworkopt);
        return;
    }
    //
    //     Quick return if possible
    //
    if (min(m, n) == 0) {
        work[1 - 1] = COMPLEX(lworkopt);
        return;
    }
    //
    //     (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in
    //     of M-by-M orthogonal matrix Q_in, which is implicitly stored in
    //     the subdiagonal part of input array A and in the input array T.
    //     Perform by the following operation using the routine Clamtsqr.
    //
    //         Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix,
    //                        ( 0 )        0 is a (M-N)-by-N zero matrix.
    //
    //     (1a) Form M-by-N matrix in the array WORK(1:LDC*N) with ones
    //     on the diagonal and zeros elsewhere.
    //
    const COMPLEX czero = (0.0, 0.0);
    const COMPLEX cone = (1.0, 0.0);
    Claset("F", m, n, czero, cone, work, ldc);
    //
    //     (1b)  On input, WORK(1:LDC*N) stores ( I );
    //                                          ( 0 )
    //
    //           On output, WORK(1:LDC*N) stores Q1_in.
    //
    INTEGER iinfo = 0;
    Clamtsqr("L", "N", m, n, n, mb, nblocal, a, lda, t, ldt, work, ldc, work[(lc + 1) - 1], lw, iinfo);
    //
    //     (2) Copy the result from the part of the work array (1:M,1:N)
    //     the output array A(1:M,1:N) column-by-column.
    //
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Ccopy(m, work[((j - 1) * ldc + 1) - 1], 1, a[(j - 1) * lda], 1);
    }
    //
    work[1 - 1] = COMPLEX(lworkopt);
    //
    //     End of Cungtsqr
    //
}

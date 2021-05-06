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

void Cgetsqrhrt(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER const nb1, INTEGER const nb2, COMPLEX *a, INTEGER const lda, COMPLEX *t, INTEGER const ldt, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    //     Test the input arguments
    //
    info = 0;
    bool lquery = lwork == -1;
    INTEGER nb1local = 0;
    INTEGER num_all_row_blocks = 0;
    INTEGER lwt = 0;
    INTEGER ldwt = 0;
    INTEGER lw1 = 0;
    INTEGER lw2 = 0;
    INTEGER lworkopt = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0 || m < n) {
        info = -2;
    } else if (mb1 <= n) {
        info = -3;
    } else if (nb1 < 1) {
        info = -4;
    } else if (nb2 < 1) {
        info = -5;
    } else if (lda < max((INTEGER)1, m)) {
        info = -7;
    } else if (ldt < max({(INTEGER)1, min(nb2, n)})) {
        info = -9;
    } else {
        //
        //        This workspace is used to store array:
        //        a) Matrix T and WORK for Clatsqr;
        //        b) N-by-N upper-triangular factor R_tsqr;
        //        c) Matrix T and array WORK for Cungtsqr_row;
        //        d) Diagonal D for Cunhr_col.
        //
        if (lwork < n * n + 1 && !lquery) {
            info = -11;
        } else {
            //
            //           Set block size for column blocks
            //
            nb1local = min(nb1, n);
            //
            num_all_row_blocks = max((INTEGER)1, ceil(castREAL(m - n) / castREAL(mb1 - n)));
            //
            //           T array in TSQR.
            //
            lwt = num_all_row_blocks * n * nb1local;
            //
            ldwt = nb1local;
            //
            //           Length of TSQR work array
            //
            lw1 = nb1local * n;
            //
            //           Length of Cungtsqr_row work array.
            //
            lw2 = nb1local * max(nb1local, (n - nb1local));
            //
            lworkopt = max({lwt + lw1, max(lwt + n * n + lw2, lwt + n * n + n)});
            //
            if ((lwork < max((INTEGER)1, lworkopt)) && (!lquery)) {
                info = -11;
            }
            //
        }
    }
    //
    //     Handle error in the input parameters and return workspace query.
    //
    if (info != 0) {
        Mxerbla("Cgetsqrhrt", -info);
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
    INTEGER nb2local = min(nb2, n);
    //
    //     (1) Perform TSQR-factorization of the M-by-N matrix A.
    //
    INTEGER iinfo = 0;
    Clatsqr(m, n, mb1, nb1local, a, lda, work, ldwt, &work[(lwt + 1) - 1], lw1, iinfo);
    //
    //     (2) Copy the factor R_tsqr stored in the upper-triangular part
    //         of A into the square matrix in the work array
    //         WORK(LWT+1:LWT+N*N) column-by-column.
    //
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Ccopy(j, &a[(j - 1) * lda], 1, &work[(lwt + n * (j - 1) + 1) - 1], 1);
    }
    //
    //     (3) Generate a M-by-N matrix Q with orthonormal columns from
    //     the result stored below the diagonal in the array A in place.
    //
    Cungtsqr_row(m, n, mb1, nb1local, a, lda, work, ldwt, &work[(lwt + n * n + 1) - 1], lw2, iinfo);
    //
    //     (4) Perform the reconstruction of Householder vectors from
    //     the matrix Q (stored in A) in place.
    //
    Cunhr_col(m, n, nb2local, a, lda, t, ldt, &work[(lwt + n * n + 1) - 1], iinfo);
    //
    //     (5) Copy the factor R_tsqr stored in the square matrix in the
    //     work array WORK(LWT+1:LWT+N*N) into the upper-triangular
    //     part of A.
    //
    //     (6) Compute from R_tsqr the factor R_hr corresponding to
    //     the reconstructed Householder vectors, i.e. R_hr = S * R_tsqr.
    //     This multiplication by the sign matrix S on the left means
    //     changing the sign of I-th row of the matrix R_tsqr according
    //     to sign of the I-th diagonal element DIAG(I) of the matrix S.
    //     DIAG is stored in WORK( LWT+N*N+1 ) from the Cunhr_col output.
    //
    //     (5) and (6) can be combined in a single loop, so the rows in A
    //     are accessed only once.
    //
    INTEGER i = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    for (i = 1; i <= n; i = i + 1) {
        if (work[(lwt + n * n + i) - 1] == -cone) {
            for (j = i; j <= n; j = j + 1) {
                a[(i - 1) + (j - 1) * lda] = -cone * work[(lwt + n * (j - 1) + i) - 1];
            }
        } else {
            Ccopy(n - i + 1, &work[(lwt + n * (i - 1) + i) - 1], n, &a[(i - 1) + (i - 1) * lda], lda);
        }
    }
    //
    work[1 - 1] = COMPLEX(lworkopt);
    //
    //     End of Cgetsqrhrt
    //
}

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

void Claunhr_col_getrfnp(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *d, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Claunhr_col_getrfnp", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (min(m, n) == 0) {
        return;
    }
    //
    //     Determine the block size for this environment.
    //
    INTEGER nb = iMlaenv(1, "Claunhr_col_getrfnp", " ", m, n, -1, -1);
    //
    INTEGER j = 0;
    INTEGER jb = 0;
    INTEGER iinfo = 0;
    const COMPLEX cone = (1.0, 0.0);
    if (nb <= 1 || nb >= min(m, n)) {
        //
        //        Use unblocked code.
        //
        Claunhr_col_getrfnp2(m, n, a, lda, d, info);
    } else {
        //
        //        Use blocked code.
        //
        for (j = 1; j <= min(m, n); j = j + nb) {
            jb = min(min(m, n) - j + 1, nb);
            //
            //           Factor diagonal and subdiagonal blocks.
            //
            Claunhr_col_getrfnp2(m - j + 1, jb, &a[(j - 1) + (j - 1) * lda], lda, &d[j - 1], iinfo);
            //
            if (j + jb <= n) {
                //
                //              Compute block row of U.
                //
                Ctrsm("Left", "Lower", "No transpose", "Unit", jb, n - j - jb + 1, cone, &a[(j - 1) + (j - 1) * lda], lda, &a[(j - 1) + ((j + jb) - 1) * lda], lda);
                if (j + jb <= m) {
                    //
                    //                 Update trailing submatrix.
                    //
                    Cgemm("No transpose", "No transpose", m - j - jb + 1, n - j - jb + 1, jb, -cone, &a[((j + jb) - 1) + (j - 1) * lda], lda, &a[(j - 1) + ((j + jb) - 1) * lda], lda, cone, &a[((j + jb) - 1) + ((j + jb) - 1) * lda], lda);
                }
            }
        }
    }
    //
    //     End of Claunhr_col_getrfnp
    //
}

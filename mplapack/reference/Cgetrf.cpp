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

void Cgetrf(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info) {
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
        Mxerbla("Cgetrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Determine the block size for this environment.
    //
    INTEGER nb = iMlaenv(1, "Cgetrf", " ", m, n, -1, -1);
    INTEGER j = 0;
    INTEGER jb = 0;
    INTEGER iinfo = 0;
    INTEGER i = 0;
    const COMPLEX one = (1.0, 0.0);
    if (nb <= 1 || nb >= min(m, n)) {
        //
        //        Use unblocked code.
        //
        Cgetrf2(m, n, a, lda, ipiv, info);
    } else {
        //
        //        Use blocked code.
        //
        for (j = 1; j <= min(m, n); j = j + nb) {
            jb = min(min(m, n) - j + 1, nb);
            //
            //           Factor diagonal and subdiagonal blocks and test for exact
            //           singularity.
            //
            Cgetrf2(m - j + 1, jb, &a[(j - 1) + (j - 1) * lda], lda, &ipiv[j - 1], iinfo);
            //
            //           Adjust INFO and the pivot indices.
            //
            if (info == 0 && iinfo > 0) {
                info = iinfo + j - 1;
            }
            for (i = j; i <= min(m, j + jb - 1); i = i + 1) {
                ipiv[i - 1] += j - 1;
            }
            //
            //           Apply interchanges to columns 1:J-1.
            //
            Claswp(j - 1, a, lda, j, j + jb - 1, ipiv, 1);
            //
            if (j + jb <= n) {
                //
                //              Apply interchanges to columns J+JB:N.
                //
                Claswp(n - j - jb + 1, &a[((j + jb) - 1) * lda], lda, j, j + jb - 1, ipiv, 1);
                //
                //              Compute block row of U.
                //
                Ctrsm("Left", "Lower", "No transpose", "Unit", jb, n - j - jb + 1, one, &a[(j - 1) + (j - 1) * lda], lda, &a[(j - 1) + ((j + jb) - 1) * lda], lda);
                if (j + jb <= m) {
                    //
                    //                 Update trailing submatrix.
                    //
                    Cgemm("No transpose", "No transpose", m - j - jb + 1, n - j - jb + 1, jb, -one, &a[((j + jb) - 1) + (j - 1) * lda], lda, &a[(j - 1) + ((j + jb) - 1) * lda], lda, one, &a[((j + jb) - 1) + ((j + jb) - 1) * lda], lda);
                }
            }
        }
    }
    //
    //     End of Cgetrf
    //
}

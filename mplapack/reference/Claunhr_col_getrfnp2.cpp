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

REAL abs1(COMPLEX ff) { return max(abs(ff.real()), abs(ff.imag())); }

void Claunhr_col_getrfnp2(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *d, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX z = 0.0;
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
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
        Mxerbla("Claunhr_col_getrfnp2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (min(m, n) == 0) {
        return;
    }
    //
    const REAL one = 1.0;
    REAL sfmin = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER i = 0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER iinfo = 0;
    if (m == 1) {
        //
        //        One row case, (also recursion termination case),
        //        use unblocked code
        //
        //        Transfer the sign
        //
        d[1 - 1] = COMPLEX(-sign(one, a[(1 - 1)].real()), 0.0);
        //
        //        Construct the row of U
        //
        a[(1 - 1)] = a[(1 - 1)] - d[1 - 1];
        //
    } else if (n == 1) {
        //
        //        One column case, (also recursion termination case),
        //        use unblocked code
        //
        //        Transfer the sign
        //
        d[1 - 1] = COMPLEX(-sign(one, a[(1 - 1)].real()), 0.0);
        //
        //        Construct the row of U
        //
        a[(1 - 1)] = a[(1 - 1)] - d[1 - 1];
        //
        //        Scale the elements 2:M of the column
        //
        //        Determine machine safe minimum
        //
        sfmin = Rlamch("S");
        //
        //        Construct the subdiagonal elements of L
        //
        if (abs1(a[(1 - 1)]) >= sfmin) {
            Cscal(m - 1, cone / a[(1 - 1)], &a[(2 - 1)], 1);
        } else {
            for (i = 2; i <= m; i = i + 1) {
                a[(i - 1)] = a[(i - 1)] / a[(1 - 1)];
            }
        }
        //
    } else {
        //
        //        Divide the matrix B into four submatrices
        //
        n1 = min(m, n) / 2;
        n2 = n - n1;
        //
        //        Factor B11, recursive call
        //
        Claunhr_col_getrfnp2(n1, n1, a, lda, d, iinfo);
        //
        //        Solve for B21
        //
        Ctrsm("R", "U", "N", "N", m - n1, n1, cone, a, lda, &a[((n1 + 1) - 1)], lda);
        //
        //        Solve for B12
        //
        Ctrsm("L", "L", "N", "U", n1, n2, cone, a, lda, &a[((n1 + 1) - 1) * lda], lda);
        //
        //        Update B22, i.e. compute the Schur complement
        //        B22 := B22 - B21*B12
        //
        Cgemm("N", "N", m - n1, n2, n1, -cone, &a[((n1 + 1) - 1)], lda, &a[((n1 + 1) - 1) * lda], lda, cone, &a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda);
        //
        //        Factor B22, recursive call
        //
        Claunhr_col_getrfnp2(m - n1, n2, &a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda, &d[(n1 + 1) - 1], iinfo);
        //
    }
    //
    //     End of Claunhr_col_getrfnp2
    //
}

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

void Cunhr_col(INTEGER const &m, INTEGER const &n, INTEGER const &nb, COMPLEX *a, INTEGER const &lda, COMPLEX *t, INTEGER const &ldt, COMPLEX *d, INTEGER &info) {
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
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0 || n > m) {
        info = -2;
    } else if (nb < 1) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldt < max((INTEGER)1, min(nb, n))) {
        info = -7;
    }
    //
    //     Handle error in the input parameters.
    //
    if (info != 0) {
        Mxerbla("Cunhr_col", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (min(m, n) == 0) {
        return;
    }
    //
    //     On input, the M-by-N matrix A contains the unitary
    //     M-by-N matrix Q_in.
    //
    //     (1) Compute the unit lower-trapezoidal V (ones on the diagonal
    //     are not stored) by performing the "modified" LU-decomposition.
    //
    //     Q_in - ( S ) = V * U = ( V1 ) * U,
    //            ( 0 )           ( V2 )
    //
    //     where 0 is an (M-N)-by-N zero matrix.
    //
    //     (1-1) Factor V1 and U.
    //
    INTEGER iinfo = 0;
    Claunhr_col_getrfnp(n, n, a, lda, d, iinfo);
    //
    //     (1-2) Solve for V2.
    //
    const COMPLEX cone = (1.0, 0.0);
    if (m > n) {
        Ctrsm("R", "U", "N", "N", m - n, n, cone, a, lda, a[((n + 1) - 1)], lda);
    }
    //
    //     (2) Reconstruct the block reflector T stored in T(1:NB, 1:N)
    //     as a sequence of upper-triangular blocks with NB-size column
    //     blocking.
    //
    //     Loop over the column blocks of size NB of the array A(1:M,1:N)
    //     and the array T(1:NB,1:N), JB is the column index of a column
    //     block, JNB is the column block size at each step JB.
    //
    INTEGER nplusone = n + 1;
    INTEGER jb = 0;
    INTEGER jnb = 0;
    INTEGER jbtemp1 = 0;
    INTEGER j = 0;
    INTEGER jbtemp2 = 0;
    INTEGER i = 0;
    const COMPLEX czero = (0.0, 0.0);
    for (jb = 1; jb <= n; jb = jb + nb) {
        //
        //        (2-0) Determine the column block size JNB.
        //
        jnb = min(nplusone - jb, nb);
        //
        //        (2-1) Copy the upper-triangular part of the current JNB-by-JNB
        //        diagonal block U(JB) (of the N-by-N matrix U) stored
        //        in A(JB:JB+JNB-1,JB:JB+JNB-1) INTEGERo the upper-triangular part
        //        of the current JNB-by-JNB block T(1:JNB,JB:JB+JNB-1)
        //        column-by-column, total JNB*(JNB+1)/2 elements.
        //
        jbtemp1 = jb - 1;
        for (j = jb; j <= jb + jnb - 1; j = j + 1) {
            Ccopy(j - jbtemp1, a[(jb - 1) + (j - 1) * lda], 1, t[(j - 1) * ldt], 1);
        }
        //
        //        (2-2) Perform on the upper-triangular part of the current
        //        JNB-by-JNB diagonal block U(JB) (of the N-by-N matrix U) stored
        //        in T(1:JNB,JB:JB+JNB-1) the following operation in place:
        //        (-1)*U(JB)*S(JB), i.e the result will be stored in the upper-
        //        triangular part of T(1:JNB,JB:JB+JNB-1). This multiplication
        //        of the JNB-by-JNB diagonal block U(JB) by the JNB-by-JNB
        //        diagonal block S(JB) of the N-by-N sign matrix S from the
        //        right means changing the sign of each J-th column of the block
        //        U(JB) according to the sign of the diagonal element of the block
        //        S(JB), i.e. S(J,J) that is stored in the array element D(J).
        //
        for (j = jb; j <= jb + jnb - 1; j = j + 1) {
            if (d[j - 1] == cone) {
                Cscal(j - jbtemp1, -cone, t[(j - 1) * ldt], 1);
            }
        }
        //
        //        (2-3) Perform the triangular solve for the current block
        //        matrix X(JB):
        //
        //               X(JB) * (A(JB)**T) = B(JB), where:
        //
        //               A(JB)**T  is a JNB-by-JNB unit upper-triangular
        //                         coefficient block, and A(JB)=V1(JB), which
        //                         is a JNB-by-JNB unit lower-triangular block
        //                         stored in A(JB:JB+JNB-1,JB:JB+JNB-1).
        //                         The N-by-N matrix V1 is the upper part
        //                         of the M-by-N lower-trapezoidal matrix V
        //                         stored in A(1:M,1:N);
        //
        //               B(JB)     is a JNB-by-JNB  upper-triangular right-hand
        //                         side block, B(JB) = (-1)*U(JB)*S(JB), and
        //                         B(JB) is stored in T(1:JNB,JB:JB+JNB-1);
        //
        //               X(JB)     is a JNB-by-JNB upper-triangular solution
        //                         block, X(JB) is the upper-triangular block
        //                         reflector T(JB), and X(JB) is stored
        //                         in T(1:JNB,JB:JB+JNB-1).
        //
        //             In other words, we perform the triangular solve for the
        //             upper-triangular block T(JB):
        //
        //               T(JB) * (V1(JB)**T) = (-1)*U(JB)*S(JB).
        //
        //             Even though the blocks X(JB) and B(JB) are upper-
        //             triangular, the routine Ctrsm will access all JNB**2
        //             elements of the square T(1:JNB,JB:JB+JNB-1). Therefore,
        //             we need to set to zero the elements of the block
        //             T(1:JNB,JB:JB+JNB-1) below the diagonal before the call
        //             to Ctrsm.
        //
        //        (2-3a) Set the elements to zero.
        //
        jbtemp2 = jb - 2;
        for (j = jb; j <= jb + jnb - 2; j = j + 1) {
            for (i = j - jbtemp2; i <= nb; i = i + 1) {
                t[(i - 1) + (j - 1) * ldt] = czero;
            }
        }
        //
        //        (2-3b) Perform the triangular solve.
        //
        Ctrsm("R", "L", "C", "U", jnb, jnb, cone, a[(jb - 1) + (jb - 1) * lda], lda, t[(jb - 1) * ldt], ldt);
        //
    }
    //
    //     End of Cunhr_col
    //
}

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

void Cpftri(const char *transr, const char *uplo, INTEGER const &n, COMPLEX *a, INTEGER &info) {
    a(dim1(0, star));
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
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
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    bool normaltransr = Mlsame(transr, "N");
    bool lower = Mlsame(uplo, "L");
    if (!normaltransr && !Mlsame(transr, "C")) {
        info = -1;
    } else if (!lower && !Mlsame(uplo, "U")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    }
    if (info != 0) {
        Mxerbla("Cpftri", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Invert the triangular Cholesky factor U or L.
    //
    Ctftri(transr, uplo, "N", n, a, info);
    if (info > 0) {
        return;
    }
    //
    //     If N is odd, set NISODD = .TRUE.
    //     If N is even, set K = N/2 and NISODD = .FALSE.
    //
    INTEGER k = 0;
    bool nisodd = false;
    if (mod(n, 2) == 0) {
        k = n / 2;
        nisodd = false;
    } else {
        nisodd = true;
    }
    //
    //     Set N1 and N2 depending on LOWER
    //
    INTEGER n2 = 0;
    INTEGER n1 = 0;
    if (lower) {
        n2 = n / 2;
        n1 = n - n2;
    } else {
        n1 = n / 2;
        n2 = n - n1;
    }
    //
    //     Start execution of triangular matrix multiply: inv(U)*inv(U)^C or
    //     inv(L)^C*inv(L). There are eight cases.
    //
    const REAL one = 1.0;
    const COMPLEX cone = (1.0, 0.0);
    if (nisodd) {
        //
        //        N is odd
        //
        if (normaltransr) {
            //
            //           N is odd and TRANSR = 'N'
            //
            if (lower) {
                //
                //              SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:N1-1) )
                //              T1 -> a(0,0), T2 -> a(0,1), S -> a(N1,0)
                //              T1 -> a(0), T2 -> a(n), S -> a(N1)
                //
                Clauum("L", n1, a[0 - 1], n, info);
                Cherk("L", "C", n1, n2, one, a[n1 - 1], n, one, a[0 - 1], n);
                Ctrmm("L", "U", "N", "N", n2, n1, cone, a[n - 1], n, a[n1 - 1], n);
                Clauum("U", n2, a[n - 1], n, info);
                //
            } else {
                //
                //              SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:N2-1)
                //              T1 -> a(N1+1,0), T2 -> a(N1,0), S -> a(0,0)
                //              T1 -> a(N2), T2 -> a(N1), S -> a(0)
                //
                Clauum("L", n1, a[n2 - 1], n, info);
                Cherk("L", "N", n1, n2, one, a[0 - 1], n, one, a[n2 - 1], n);
                Ctrmm("R", "U", "C", "N", n1, n2, cone, a[n1 - 1], n, a[0 - 1], n);
                Clauum("U", n2, a[n1 - 1], n, info);
                //
            }
            //
        } else {
            //
            //           N is odd and TRANSR = 'C'
            //
            if (lower) {
                //
                //              SRPA for LOWER, TRANSPOSE, and N is odd
                //              T1 -> a(0), T2 -> a(1), S -> a(0+N1*N1)
                //
                Clauum("U", n1, a[0 - 1], n1, info);
                Cherk("U", "N", n1, n2, one, a[(n1 * n1) - 1], n1, one, a[0 - 1], n1);
                Ctrmm("R", "L", "N", "N", n1, n2, cone, a[1 - 1], n1, a[(n1 * n1) - 1], n1);
                Clauum("L", n2, a[1 - 1], n1, info);
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE, and N is odd
                //              T1 -> a(0+N2*N2), T2 -> a(0+N1*N2), S -> a(0)
                //
                Clauum("U", n1, a[(n2 * n2) - 1], n2, info);
                Cherk("U", "C", n1, n2, one, a[0 - 1], n2, one, a[(n2 * n2) - 1], n2);
                Ctrmm("L", "L", "C", "N", n2, n1, cone, a[(n1 * n2) - 1], n2, a[0 - 1], n2);
                Clauum("L", n2, a[(n1 * n2) - 1], n2, info);
                //
            }
            //
        }
        //
    } else {
        //
        //        N is even
        //
        if (normaltransr) {
            //
            //           N is even and TRANSR = 'N'
            //
            if (lower) {
                //
                //              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
                //              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
                //              T1 -> a(1), T2 -> a(0), S -> a(k+1)
                //
                Clauum("L", k, a[1 - 1], n + 1, info);
                Cherk("L", "C", k, k, one, a[(k + 1) - 1], n + 1, one, a[1 - 1], n + 1);
                Ctrmm("L", "U", "N", "N", k, k, cone, a[0 - 1], n + 1, a[(k + 1) - 1], n + 1);
                Clauum("U", k, a[0 - 1], n + 1, info);
                //
            } else {
                //
                //              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
                //              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
                //              T1 -> a(k+1), T2 -> a(k), S -> a(0)
                //
                Clauum("L", k, a[(k + 1) - 1], n + 1, info);
                Cherk("L", "N", k, k, one, a[0 - 1], n + 1, one, a[(k + 1) - 1], n + 1);
                Ctrmm("R", "U", "C", "N", k, k, cone, a[k - 1], n + 1, a[0 - 1], n + 1);
                Clauum("U", k, a[k - 1], n + 1, info);
                //
            }
            //
        } else {
            //
            //           N is even and TRANSR = 'C'
            //
            if (lower) {
                //
                //              SRPA for LOWER, TRANSPOSE, and N is even (see paper)
                //              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1),
                //              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
                //
                Clauum("U", k, a[k - 1], k, info);
                Cherk("U", "N", k, k, one, a[(k * (k + 1)) - 1], k, one, a[k - 1], k);
                Ctrmm("R", "L", "N", "N", k, k, cone, a[0 - 1], k, a[(k * (k + 1)) - 1], k);
                Clauum("L", k, a[0 - 1], k, info);
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE, and N is even (see paper)
                //              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0),
                //              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
                //
                Clauum("U", k, a[(k * (k + 1)) - 1], k, info);
                Cherk("U", "C", k, k, one, a[0 - 1], k, one, a[(k * (k + 1)) - 1], k);
                Ctrmm("L", "L", "C", "N", k, k, cone, a[(k * k) - 1], k, a[0 - 1], k);
                Clauum("L", k, a[(k * k) - 1], k, info);
                //
            }
            //
        }
        //
    }
    //
    //     End of Cpftri
    //
}

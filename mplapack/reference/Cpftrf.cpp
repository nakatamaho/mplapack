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

void Cpftrf(const char *transr, const char *uplo, INTEGER const &n, COMPLEX *a, INTEGER &info) {
    a(dim1(0, star));
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
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
        Mxerbla("Cpftrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
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
    //     start execution: there are eight cases
    //
    const COMPLEX cone = (1.0, 0.0);
    const REAL one = 1.0;
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
                //             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
                //             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
                //             T1 -> a(0), T2 -> a(n), S -> a(n1)
                //
                Cpotrf("L", n1, a[0 - 1], n, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("R", "L", "C", "N", n2, n1, cone, a[0 - 1], n, a[n1 - 1], n);
                Cherk("U", "N", n2, n1, -one, a[n1 - 1], n, one, a[n - 1], n);
                Cpotrf("U", n2, a[n - 1], n, info);
                if (info > 0) {
                    info += n1;
                }
                //
            } else {
                //
                //             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
                //             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
                //             T1 -> a(n2), T2 -> a(n1), S -> a(0)
                //
                Cpotrf("L", n1, a[n2 - 1], n, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("L", "L", "N", "N", n1, n2, cone, a[n2 - 1], n, a[0 - 1], n);
                Cherk("U", "C", n2, n1, -one, a[0 - 1], n, one, a[n1 - 1], n);
                Cpotrf("U", n2, a[n1 - 1], n, info);
                if (info > 0) {
                    info += n1;
                }
                //
            }
            //
        } else {
            //
            //           N is odd and TRANSR = 'C'
            //
            if (lower) {
                //
                //              SRPA for LOWER, TRANSPOSE and N is odd
                //              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
                //              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1
                //
                Cpotrf("U", n1, a[0 - 1], n1, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("L", "U", "C", "N", n1, n2, cone, a[0 - 1], n1, a[(n1 * n1) - 1], n1);
                Cherk("L", "C", n2, n1, -one, a[(n1 * n1) - 1], n1, one, a[1 - 1], n1);
                Cpotrf("L", n2, a[1 - 1], n1, info);
                if (info > 0) {
                    info += n1;
                }
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is odd
                //              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
                //              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2
                //
                Cpotrf("U", n1, a[(n2 * n2) - 1], n2, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("R", "U", "N", "N", n2, n1, cone, a[(n2 * n2) - 1], n2, a[0 - 1], n2);
                Cherk("L", "N", n2, n1, -one, a[0 - 1], n2, one, a[(n1 * n2) - 1], n2);
                Cpotrf("L", n2, a[(n1 * n2) - 1], n2, info);
                if (info > 0) {
                    info += n1;
                }
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
                Cpotrf("L", k, a[1 - 1], n + 1, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("R", "L", "C", "N", k, k, cone, a[1 - 1], n + 1, a[(k + 1) - 1], n + 1);
                Cherk("U", "N", k, k, -one, a[(k + 1) - 1], n + 1, one, a[0 - 1], n + 1);
                Cpotrf("U", k, a[0 - 1], n + 1, info);
                if (info > 0) {
                    info += k;
                }
                //
            } else {
                //
                //              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
                //              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
                //              T1 -> a(k+1), T2 -> a(k), S -> a(0)
                //
                Cpotrf("L", k, a[(k + 1) - 1], n + 1, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("L", "L", "N", "N", k, k, cone, a[(k + 1) - 1], n + 1, a[0 - 1], n + 1);
                Cherk("U", "C", k, k, -one, a[0 - 1], n + 1, one, a[k - 1], n + 1);
                Cpotrf("U", k, a[k - 1], n + 1, info);
                if (info > 0) {
                    info += k;
                }
                //
            }
            //
        } else {
            //
            //           N is even and TRANSR = 'C'
            //
            if (lower) {
                //
                //              SRPA for LOWER, TRANSPOSE and N is even (see paper)
                //              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
                //              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
                //
                Cpotrf("U", k, a[(0 + k) - 1], k, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("L", "U", "C", "N", k, k, cone, a[k - 1], n1, a[(k * (k + 1)) - 1], k);
                Cherk("L", "C", k, k, -one, a[(k * (k + 1)) - 1], k, one, a[0 - 1], k);
                Cpotrf("L", k, a[0 - 1], k, info);
                if (info > 0) {
                    info += k;
                }
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is even (see paper)
                //              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
                //              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
                //
                Cpotrf("U", k, a[(k * (k + 1)) - 1], k, info);
                if (info > 0) {
                    return;
                }
                Ctrsm("R", "U", "N", "N", k, k, cone, a[(k * (k + 1)) - 1], k, a[0 - 1], k);
                Cherk("L", "N", k, k, -one, a[0 - 1], k, one, a[(k * k) - 1], k);
                Cpotrf("L", k, a[(k * k) - 1], k, info);
                if (info > 0) {
                    info += k;
                }
                //
            }
            //
        }
        //
    }
    //
    //     End of Cpftrf
    //
}

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

void Rtftri(const char *transr, const char *uplo, const char *diag, INTEGER const &n, REAL *a, INTEGER &info) {
    a(dim1(0, star));
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
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    bool normaltransr = Mlsame(transr, "N");
    bool lower = Mlsame(uplo, "L");
    if (!normaltransr && !Mlsame(transr, "T")) {
        info = -1;
    } else if (!lower && !Mlsame(uplo, "U")) {
        info = -2;
    } else if (!Mlsame(diag, "N") && !Mlsame(diag, "U")) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Rtftri", -info);
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
                Rtrtri("L", diag, n1, a[0 - 1], n, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "L", "N", diag, n2, n1, -one, a[0 - 1], n, a[n1 - 1], n);
                Rtrtri("U", diag, n2, a[n - 1], n, info);
                if (info > 0) {
                    info += n1;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "U", "T", diag, n2, n1, one, a[n - 1], n, a[n1 - 1], n);
                //
            } else {
                //
                //             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
                //             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
                //             T1 -> a(n2), T2 -> a(n1), S -> a(0)
                //
                Rtrtri("L", diag, n1, a[n2 - 1], n, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "L", "T", diag, n1, n2, -one, a[n2 - 1], n, a[0 - 1], n);
                Rtrtri("U", diag, n2, a[n1 - 1], n, info);
                if (info > 0) {
                    info += n1;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "U", "N", diag, n1, n2, one, a[n1 - 1], n, a[0 - 1], n);
                //
            }
            //
        } else {
            //
            //           N is odd and TRANSR = 'T'
            //
            if (lower) {
                //
                //              SRPA for LOWER, TRANSPOSE and N is odd
                //              T1 -> a(0), T2 -> a(1), S -> a(0+n1*n1)
                //
                Rtrtri("U", diag, n1, a[0 - 1], n1, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "U", "N", diag, n1, n2, -one, a[0 - 1], n1, a[(n1 * n1) - 1], n1);
                Rtrtri("L", diag, n2, a[1 - 1], n1, info);
                if (info > 0) {
                    info += n1;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "L", "T", diag, n1, n2, one, a[1 - 1], n1, a[(n1 * n1) - 1], n1);
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is odd
                //              T1 -> a(0+n2*n2), T2 -> a(0+n1*n2), S -> a(0)
                //
                Rtrtri("U", diag, n1, a[(n2 * n2) - 1], n2, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "U", "T", diag, n2, n1, -one, a[(n2 * n2) - 1], n2, a[0 - 1], n2);
                Rtrtri("L", diag, n2, a[(n1 * n2) - 1], n2, info);
                if (info > 0) {
                    info += n1;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "L", "N", diag, n2, n1, one, a[(n1 * n2) - 1], n2, a[0 - 1], n2);
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
                Rtrtri("L", diag, k, a[1 - 1], n + 1, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "L", "N", diag, k, k, -one, a[1 - 1], n + 1, a[(k + 1) - 1], n + 1);
                Rtrtri("U", diag, k, a[0 - 1], n + 1, info);
                if (info > 0) {
                    info += k;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "U", "T", diag, k, k, one, a[0 - 1], n + 1, a[(k + 1) - 1], n + 1);
                //
            } else {
                //
                //              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
                //              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
                //              T1 -> a(k+1), T2 -> a(k), S -> a(0)
                //
                Rtrtri("L", diag, k, a[(k + 1) - 1], n + 1, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "L", "T", diag, k, k, -one, a[(k + 1) - 1], n + 1, a[0 - 1], n + 1);
                Rtrtri("U", diag, k, a[k - 1], n + 1, info);
                if (info > 0) {
                    info += k;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "U", "N", diag, k, k, one, a[k - 1], n + 1, a[0 - 1], n + 1);
            }
        } else {
            //
            //           N is even and TRANSR = 'T'
            //
            if (lower) {
                //
                //              SRPA for LOWER, TRANSPOSE and N is even (see paper)
                //              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
                //              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
                //
                Rtrtri("U", diag, k, a[k - 1], k, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "U", "N", diag, k, k, -one, a[k - 1], k, a[(k * (k + 1)) - 1], k);
                Rtrtri("L", diag, k, a[0 - 1], k, info);
                if (info > 0) {
                    info += k;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "L", "T", diag, k, k, one, a[0 - 1], k, a[(k * (k + 1)) - 1], k);
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is even (see paper)
                //              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
                //              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
                //
                Rtrtri("U", diag, k, a[(k * (k + 1)) - 1], k, info);
                if (info > 0) {
                    return;
                }
                Rtrmm("R", "U", "T", diag, k, k, -one, a[(k * (k + 1)) - 1], k, a[0 - 1], k);
                Rtrtri("L", diag, k, a[(k * k) - 1], k, info);
                if (info > 0) {
                    info += k;
                }
                if (info > 0) {
                    return;
                }
                Rtrmm("L", "L", "N", diag, k, k, one, a[(k * k) - 1], k, a[0 - 1], k);
            }
        }
    }
    //
    //     End of Rtftri
    //
}

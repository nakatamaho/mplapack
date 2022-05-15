/*
 * Copyright (c) 2021-2022
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

void Rpftrf(const char *transr, const char *uplo, INTEGER const n, REAL *a, INTEGER &info) {
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
    } else if (n < 0) {
        info = -3;
    }
    if (info != 0) {
        Mxerbla("Rpftrf", -info);
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
                Rpotrf("L", n1, &a[0], n, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("R", "L", "T", "N", n2, n1, one, &a[0], n, &a[n1], n);
                Rsyrk("U", "N", n2, n1, -one, &a[n1], n, one, &a[n], n);
                Rpotrf("U", n2, &a[n], n, info);
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
                Rpotrf("L", n1, &a[n2], n, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("L", "L", "N", "N", n1, n2, one, &a[n2], n, &a[0], n);
                Rsyrk("U", "T", n2, n1, -one, &a[0], n, one, &a[n1], n);
                Rpotrf("U", n2, &a[n1], n, info);
                if (info > 0) {
                    info += n1;
                }
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
                //              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
                //              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1
                //
                Rpotrf("U", n1, &a[0], n1, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("L", "U", "T", "N", n1, n2, one, &a[0], n1, &a[(n1 * n1)], n1);
                Rsyrk("L", "T", n2, n1, -one, &a[(n1 * n1)], n1, one, &a[1], n1);
                Rpotrf("L", n2, &a[1], n1, info);
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
                Rpotrf("U", n1, &a[(n2 * n2)], n2, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("R", "U", "N", "N", n2, n1, one, &a[(n2 * n2)], n2, &a[0 - 1], n2);
                Rsyrk("L", "N", n2, n1, -one, &a[0], n2, one, &a[(n1 * n2)], n2);
                Rpotrf("L", n2, &a[(n1 * n2)], n2, info);
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
                Rpotrf("L", k, &a[1], n + 1, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("R", "L", "T", "N", k, k, one, &a[1 ], n + 1, &a[(k + 1) ], n + 1);
                Rsyrk("U", "N", k, k, -one, &a[(k + 1)], n + 1, one, &a[0], n + 1);
                Rpotrf("U", k, &a[0], n + 1, info);
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
                Rpotrf("L", k, &a[(k + 1)], n + 1, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("L", "L", "N", "N", k, k, one, &a[(k + 1)], n + 1, &a[0], n + 1);
                Rsyrk("U", "T", k, k, -one, &a[0], n + 1, one, &a[k], n + 1);
                Rpotrf("U", k, &a[k], n + 1, info);
                if (info > 0) {
                    info += k;
                }
                //
            }
            //
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
                Rpotrf("U", k, &a[(0 + k)], k, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("L", "U", "T", "N", k, k, one, &a[k], n1, &a[(k * (k + 1))], k);
                Rsyrk("L", "T", k, k, -one, &a[(k * (k + 1))], k, one, &a[0], k);
                Rpotrf("L", k, &a[0], k, info);
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
                Rpotrf("U", k, &a[(k * (k + 1))], k, info);
                if (info > 0) {
                    return;
                }
                Rtrsm("R", "U", "N", "N", k, k, one, &a[(k * (k + 1))], k, &a[0], k);
                Rsyrk("L", "N", k, k, -one, &a[0], k, one, &a[(k * k)], k);
                Rpotrf("L", k, &a[(k * k)], k, info);
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
    //     End of Rpftrf
    //
}

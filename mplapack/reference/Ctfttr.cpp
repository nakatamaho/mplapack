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

void Ctfttr(const char *transr, const char *uplo, INTEGER const n, COMPLEX *arf, COMPLEX *a, INTEGER const lda, INTEGER &info) {
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
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Ctfttr", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        if (n == 1) {
            if (normaltransr) {
                a[0] = arf[0];
            } else {
                a[0] = conj(arf[0]);
            }
        }
        return;
    }
    //
    //     Size of array ARF(1:2,0:nt-1)
    //
    INTEGER nt = n * (n + 1) / 2;
    //
    //     set N1 and N2 depending on LOWER: for N even N1=N2=K
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
    //     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2.
    //     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is
    //     N--by--(N+1)/2.
    //
    INTEGER k = 0;
    bool nisodd = false;
    INTEGER np1x2 = 0;
    INTEGER nx2 = 0;
    if (mod(n, 2) == 0) {
        k = n / 2;
        nisodd = false;
        if (!lower) {
            np1x2 = n + n + 2;
        }
    } else {
        nisodd = true;
        if (!lower) {
            nx2 = n + n;
        }
    }
    //
    INTEGER ij = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER l = 0;
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
                //             T1 -> a(0), T2 -> a(n), S -> a(n1); lda=n
                //
                ij = 0;
                for (j = 0; j <= n2; j = j + 1) {
                    for (i = n1; i <= n2 + j; i = i + 1) {
                        a[(n2 + j) + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                    for (i = j; i <= n - 1; i = i + 1) {
                        a[i + j * lda] = arf[ij];
                        ij++;
                    }
                }
                //
            } else {
                //
                //             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
                //             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
                //             T1 -> a(n2), T2 -> a(n1), S -> a(0); lda=n
                //
                ij = nt - n;
                for (j = n - 1; j >= n1; j = j - 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        a[i + j * lda] = arf[ij];
                        ij++;
                    }
                    for (l = j - n1; l <= n1 - 1; l = l + 1) {
                        a[(j - n1) + l * lda] = conj(arf[ij]);
                        ij++;
                    }
                    ij = ij - nx2;
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
                //              T1 -> A(0+0) , T2 -> A(1+0) , S -> A(0+n1*n1); lda=n1
                //
                ij = 0;
                for (j = 0; j <= n2 - 1; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        a[j + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                    for (i = n1 + j; i <= n - 1; i = i + 1) {
                        a[i + (n1 + j) * lda] = arf[ij];
                        ij++;
                    }
                }
                for (j = n2; j <= n - 1; j = j + 1) {
                    for (i = 0; i <= n1 - 1; i = i + 1) {
                        a[j + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                }
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is odd
                //              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
                //              T1 -> A(n2*n2), T2 -> A(n1*n2), S -> A(0); lda = n2
                //
                ij = 0;
                for (j = 0; j <= n1; j = j + 1) {
                    for (i = n1; i <= n - 1; i = i + 1) {
                        a[j + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                }
                for (j = 0; j <= n1 - 1; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        a[i + j * lda] = arf[ij];
                        ij++;
                    }
                    for (l = n2 + j; l <= n - 1; l = l + 1) {
                        a[(n2 + j) + l * lda] = conj(arf[ij]);
                        ij++;
                    }
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
                //              T1 -> a(1), T2 -> a(0), S -> a(k+1); lda=n+1
                //
                ij = 0;
                for (j = 0; j <= k - 1; j = j + 1) {
                    for (i = k; i <= k + j; i = i + 1) {
                        a[(k + j) + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                    for (i = j; i <= n - 1; i = i + 1) {
                        a[i + j * lda] = arf[ij];
                        ij++;
                    }
                }
                //
            } else {
                //
                //              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
                //              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
                //              T1 -> a(k+1), T2 -> a(k), S -> a(0); lda=n+1
                //
                ij = nt - n - 1;
                for (j = n - 1; j >= k; j = j - 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        a[i + j * lda] = arf[ij];
                        ij++;
                    }
                    for (l = j - k; l <= k - 1; l = l + 1) {
                        a[(j - k) + l * lda] = conj(arf[ij]);
                        ij++;
                    }
                    ij = ij - np1x2;
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
                //              SRPA for LOWER, TRANSPOSE and N is even (see paper, A=B)
                //              T1 -> A(0,1) , T2 -> A(0,0) , S -> A(0,k+1) :
                //              T1 -> A(0+k) , T2 -> A(0+0) , S -> A(0+k*(k+1)); lda=k
                //
                ij = 0;
                j = k;
                for (i = k; i <= n - 1; i = i + 1) {
                    a[i + j * lda] = arf[ij];
                    ij++;
                }
                for (j = 0; j <= k - 2; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        a[j + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                    for (i = k + 1 + j; i <= n - 1; i = i + 1) {
                        a[i + (k + 1 + j) * lda] = arf[ij];
                        ij++;
                    }
                }
                for (j = k - 1; j <= n - 1; j = j + 1) {
                    for (i = 0; i <= k - 1; i = i + 1) {
                        a[j + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                }
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is even (see paper, A=B)
                //              T1 -> A(0,k+1) , T2 -> A(0,k) , S -> A(0,0)
                //              T1 -> A(0+k*(k+1)) , T2 -> A(0+k*k) , S -> A(0+0)); lda=k
                //
                ij = 0;
                for (j = 0; j <= k; j = j + 1) {
                    for (i = k; i <= n - 1; i = i + 1) {
                        a[j + i * lda] = conj(arf[ij]);
                        ij++;
                    }
                }
                for (j = 0; j <= k - 2; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        a[i + j * lda] = arf[ij];
                        ij++;
                    }
                    for (l = k + 1 + j; l <= n - 1; l = l + 1) {
                        a[(k + 1 + j) + l * lda] = conj(arf[ij]);
                        ij++;
                    }
                }
                //
                //              Note that here J = K-1
                //
                for (i = 0; i <= j; i = i + 1) {
                    a[i + j * lda] = arf[ij];
                    ij++;
                }
                //
            }
            //
        }
        //
    }
    //
    //     End of Ctfttr
    //
}

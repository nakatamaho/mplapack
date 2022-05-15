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

void Rtfttp(const char *transr, const char *uplo, INTEGER const n, REAL *arf, REAL *ap, INTEGER &info) {
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
        Mxerbla("Rtfttp", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    if (n == 1) {
        if (normaltransr) {
            ap[0 - 1] = arf[0 - 1];
        } else {
            ap[0 - 1] = arf[0 - 1];
        }
        return;
    }
    //
    //     Size of array ARF(0:NT-1)
    //
    INTEGER nt = n * (n + 1) / 2;
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
    //     If N is odd, set NISODD = .TRUE.
    //     If N is even, set K = N/2 and NISODD = .FALSE.
    //
    //     set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe)
    //     where noe = 0 if n is even, noe = 1 if n is odd
    //
    INTEGER k = 0;
    bool nisodd = false;
    INTEGER lda = 0;
    if (mod(n, 2) == 0) {
        k = n / 2;
        nisodd = false;
        lda = n + 1;
    } else {
        nisodd = true;
        lda = n;
    }
    //
    //     ARF^C has lda rows and n+1-noe cols
    //
    if (!normaltransr) {
        lda = (n + 1) / 2;
    }
    //
    //     start execution: there are eight cases
    //
    INTEGER ijp = 0;
    INTEGER jp = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER ij = 0;
    INTEGER js = 0;
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
                //             T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n
                //
                ijp = 0;
                jp = 0;
                for (j = 0; j <= n2; j = j + 1) {
                    for (i = j; i <= n - 1; i = i + 1) {
                        ij = i + jp;
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    jp += lda;
                }
                for (i = 0; i <= n2 - 1; i = i + 1) {
                    for (j = 1 + i; j <= n2; j = j + 1) {
                        ij = i + j * lda;
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                }
                //
            } else {
                //
                //             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
                //             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
                //             T1 -> a(n2), T2 -> a(n1), S -> a(0)
                //
                ijp = 0;
                for (j = 0; j <= n1 - 1; j = j + 1) {
                    ij = n2 + j;
                    for (i = 0; i <= j; i = i + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                        ij += lda;
                    }
                }
                js = 0;
                for (j = n1; j <= n - 1; j = j + 1) {
                    ij = js;
                    for (ij = js; ij <= js + j; ij = ij + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    js += lda;
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
                ijp = 0;
                for (i = 0; i <= n2; i = i + 1) {
                    for (ij = i * (lda + 1); ij <= n * lda - 1; ij = ij + lda) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                }
                js = 1;
                for (j = 0; j <= n2 - 1; j = j + 1) {
                    for (ij = js; ij <= js + n2 - j - 1; ij = ij + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    js += lda + 1;
                }
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is odd
                //              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
                //              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2
                //
                ijp = 0;
                js = n2 * lda;
                for (j = 0; j <= n1 - 1; j = j + 1) {
                    for (ij = js; ij <= js + j; ij = ij + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    js += lda;
                }
                for (i = 0; i <= n1; i = i + 1) {
                    for (ij = i; ij <= i + (n1 + i) * lda; ij = ij + lda) {
                        ap[ijp] = arf[ij];
                        ijp++;
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
                //              T1 -> a(1), T2 -> a(0), S -> a(k+1)
                //
                ijp = 0;
                jp = 0;
                for (j = 0; j <= k - 1; j = j + 1) {
                    for (i = j; i <= n - 1; i = i + 1) {
                        ij = 1 + i + jp;
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    jp += lda;
                }
                for (i = 0; i <= k - 1; i = i + 1) {
                    for (j = i; j <= k - 1; j = j + 1) {
                        ij = i + j * lda;
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                }
                //
            } else {
                //
                //              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
                //              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
                //              T1 -> a(k+1), T2 -> a(k), S -> a(0)
                //
                ijp = 0;
                for (j = 0; j <= k - 1; j = j + 1) {
                    ij = k + 1 + j;
                    for (i = 0; i <= j; i = i + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                        ij += lda;
                    }
                }
                js = 0;
                for (j = k; j <= n - 1; j = j + 1) {
                    ij = js;
                    for (ij = js; ij <= js + j; ij = ij + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    js += lda;
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
                ijp = 0;
                for (i = 0; i <= k - 1; i = i + 1) {
                    for (ij = i + (i + 1) * lda; ij <= (n + 1) * lda - 1; ij = ij + lda) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                }
                js = 0;
                for (j = 0; j <= k - 1; j = j + 1) {
                    for (ij = js; ij <= js + k - j - 1; ij = ij + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    js += lda + 1;
                }
                //
            } else {
                //
                //              SRPA for UPPER, TRANSPOSE and N is even (see paper)
                //              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
                //              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
                //
                ijp = 0;
                js = (k + 1) * lda;
                for (j = 0; j <= k - 1; j = j + 1) {
                    for (ij = js; ij <= js + j; ij = ij + 1) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                    js += lda;
                }
                for (i = 0; i <= k - 1; i = i + 1) {
                    for (ij = i; ij <= i + (k + i) * lda; ij = ij + lda) {
                        ap[ijp] = arf[ij];
                        ijp++;
                    }
                }
                //
            }
            //
        }
        //
    }
    //
    //     End of Rtfttp
    //
}

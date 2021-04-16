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

void Rsytri_3x(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *e, INTEGER *ipiv, REAL *work, INTEGER const nb, INTEGER &info) {
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
    bool upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    //
    //     Quick return if possible
    //
    if (info != 0) {
        Mxerbla("Rsytri_3x", -info);
        return;
    }
    if (n == 0) {
        return;
    }
    //
    //     Workspace got Non-diag elements of D
    //
    INTEGER k = 0;
    for (k = 1; k <= n; k = k + 1) {
        work[(k - 1)] = e[k - 1];
    }
    //
    //     Check that the diagonal matrix D is nonsingular.
    //
    const REAL zero = 0.0;
    if (upper) {
        //
        //        Upper triangular storage: examine D from bottom to top
        //
        for (info = n; info >= 1; info = info - 1) {
            if (ipiv[info - 1] > 0 && a[(info - 1) + (info - 1) * lda] == zero) {
                return;
            }
        }
    } else {
        //
        //        Lower triangular storage: examine D from top to bottom.
        //
        for (info = 1; info <= n; info = info + 1) {
            if (ipiv[info - 1] > 0 && a[(info - 1) + (info - 1) * lda] == zero) {
                return;
            }
        }
    }
    //
    info = 0;
    //
    //     Splitting Workspace
    //     U01 is a block ( N, NB+1 )
    //     The first element of U01 is in WORK( 1, 1 )
    //     U11 is a block ( NB+1, NB+1 )
    //     The first element of U11 is in WORK( N+1, 1 )
    //
    INTEGER u11 = n;
    //
    //     INVD is a block ( N, 2 )
    //     The first element of INVD is in WORK( 1, INVD )
    //
    INTEGER invd = nb + 2;
    //
    const REAL one = 1.0;
    REAL t = 0.0;
    REAL ak = 0.0;
    REAL akp1 = 0.0;
    REAL akkp1 = 0.0;
    REAL d = 0.0;
    INTEGER cut = 0;
    INTEGER nnb = 0;
    INTEGER icount = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    REAL u01_i_j = 0.0;
    REAL u01_ip1_j = 0.0;
    REAL u11_i_j = 0.0;
    REAL u11_ip1_j = 0.0;
    INTEGER ip = 0;
    INTEGER ldwork = n + nb + 1;
    if (upper) {
        //
        //        Begin Upper
        //
        //        invA = P * inv(U**T) * inv(D) * inv(U) * P**T.
        //
        Rtrtri(uplo, "U", n, a, lda, info);
        //
        //        inv(D) and inv(D) * inv(U)
        //
        k = 1;
        while (k <= n) {
            if (ipiv[k - 1] > 0) {
                //              1 x 1 diagonal NNB
                work[(k - 1) + (invd - 1) * ldwork] = one / a[(k - 1) + (k - 1) * lda];
                work[(k - 1) + ((invd + 1) - 1) * ldwork] = zero;
            } else {
                //              2 x 2 diagonal NNB
                t = work[((k + 1) - 1)];
                ak = a[(k - 1) + (k - 1) * lda] / t;
                akp1 = a[((k + 1) - 1) + ((k + 1) - 1) * lda] / t;
                akkp1 = work[((k + 1) - 1)] / t;
                d = t * (ak * akp1 - one);
                work[(k - 1) + (invd - 1) * ldwork] = akp1 / d;
                work[((k + 1) - 1) + ((invd + 1) - 1) * ldwork] = ak / d;
                work[(k - 1) + ((invd + 1) - 1) * ldwork] = -akkp1 / d;
                work[((k + 1) - 1) + (invd - 1) * ldwork] = work[(k - 1) + ((invd + 1) - 1) * ldwork];
                k++;
            }
            k++;
        }
        //
        //        inv(U**T) = (inv(U))**T
        //
        //        inv(U**T) * inv(D) * inv(U)
        //
        cut = n;
        while (cut > 0) {
            nnb = nb;
            if (cut <= nnb) {
                nnb = cut;
            } else {
                icount = 0;
                //              count negative elements,
                for (i = cut + 1 - nnb; i <= cut; i = i + 1) {
                    if (ipiv[i - 1] < 0) {
                        icount++;
                    }
                }
                //              need a even number for a clear cut
                if (mod(icount, 2) == 1) {
                    nnb++;
                }
            }
            //
            cut = cut - nnb;
            //
            //           U01 Block
            //
            for (i = 1; i <= cut; i = i + 1) {
                for (j = 1; j <= nnb; j = j + 1) {
                    work[(i - 1) + (j - 1) * ldwork] = a[(i - 1) + ((cut + j) - 1) * lda];
                }
            }
            //
            //           U11 Block
            //
            for (i = 1; i <= nnb; i = i + 1) {
                work[((u11 + i) - 1) + (i - 1) * ldwork] = one;
                for (j = 1; j <= i - 1; j = j + 1) {
                    work[((u11 + i) - 1) + (j - 1) * ldwork] = zero;
                }
                for (j = i + 1; j <= nnb; j = j + 1) {
                    work[((u11 + i) - 1) + (j - 1) * ldwork] = a[((cut + i) - 1) + ((cut + j) - 1) * lda];
                }
            }
            //
            //           invD * U01
            //
            i = 1;
            while (i <= cut) {
                if (ipiv[i - 1] > 0) {
                    for (j = 1; j <= nnb; j = j + 1) {
                        work[(i - 1) + (j - 1) * ldwork] = work[(i - 1) + (invd - 1) * ldwork] * work[(i - 1) + (j - 1) * ldwork];
                    }
                } else {
                    for (j = 1; j <= nnb; j = j + 1) {
                        u01_i_j = work[(i - 1) + (j - 1) * ldwork];
                        u01_ip1_j = work[((i + 1) - 1) + (j - 1) * ldwork];
                        work[(i - 1) + (j - 1) * ldwork] = work[(i - 1) + (invd - 1) * ldwork] * u01_i_j + work[(i - 1) + ((invd + 1) - 1) * ldwork] * u01_ip1_j;
                        work[((i + 1) - 1) + (j - 1) * ldwork] = work[((i + 1) - 1) + (invd - 1) * ldwork] * u01_i_j + work[((i + 1) - 1) + ((invd + 1) - 1) * ldwork] * u01_ip1_j;
                    }
                    i++;
                }
                i++;
            }
            //
            //           invD1 * U11
            //
            i = 1;
            while (i <= nnb) {
                if (ipiv[(cut + i) - 1] > 0) {
                    for (j = i; j <= nnb; j = j + 1) {
                        work[((u11 + i) - 1) + (j - 1) * ldwork] = work[((cut + i) - 1) + (invd - 1) * ldwork] * work[((u11 + i) - 1) + (j - 1) * ldwork];
                    }
                } else {
                    for (j = i; j <= nnb; j = j + 1) {
                        u11_i_j = work[((u11 + i) - 1) + (j - 1) * ldwork];
                        u11_ip1_j = work[((u11 + i + 1) - 1) + (j - 1) * ldwork];
                        work[((u11 + i) - 1) + (j - 1) * ldwork] = work[((cut + i) - 1) + (invd - 1) * ldwork] * work[((u11 + i) - 1) + (j - 1) * ldwork] + work[((cut + i) - 1) + ((invd + 1) - 1) * ldwork] * work[((u11 + i + 1) - 1) + (j - 1) * ldwork];
                        work[((u11 + i + 1) - 1) + (j - 1) * ldwork] = work[((cut + i + 1) - 1) + (invd - 1) * ldwork] * u11_i_j + work[((cut + i + 1) - 1) + ((invd + 1) - 1) * ldwork] * u11_ip1_j;
                    }
                    i++;
                }
                i++;
            }
            //
            //           U11**T * invD1 * U11 -> U11
            //
            Rtrmm("L", "U", "T", "U", nnb, nnb, one, &a[((cut + 1) - 1) + ((cut + 1) - 1) * lda], lda, &work[((u11 + 1) - 1)], n + nb + 1);
            //
            for (i = 1; i <= nnb; i = i + 1) {
                for (j = i; j <= nnb; j = j + 1) {
                    a[((cut + i) - 1) + ((cut + j) - 1) * lda] = work[((u11 + i) - 1) + (j - 1) * ldwork];
                }
            }
            //
            //           U01**T * invD * U01 -> A( CUT+I, CUT+J )
            //
            Rgemm("T", "N", nnb, nnb, cut, one, &a[((cut + 1) - 1) * lda], lda, work, n + nb + 1, zero, &work[((u11 + 1) - 1)], n + nb + 1);
            //
            //           U11 =  U11**T * invD1 * U11 + U01**T * invD * U01
            //
            for (i = 1; i <= nnb; i = i + 1) {
                for (j = i; j <= nnb; j = j + 1) {
                    a[((cut + i) - 1) + ((cut + j) - 1) * lda] += work[((u11 + i) - 1) + (j - 1) * ldwork];
                }
            }
            //
            //           U01 =  U00**T * invD0 * U01
            //
            Rtrmm("L", uplo, "T", "U", cut, nnb, one, a, lda, work, n + nb + 1);
            //
            //           Update U01
            //
            for (i = 1; i <= cut; i = i + 1) {
                for (j = 1; j <= nnb; j = j + 1) {
                    a[(i - 1) + ((cut + j) - 1) * lda] = work[(i - 1) + (j - 1) * ldwork];
                }
            }
            //
            //           Next Block
            //
        }
        //
        //        Apply PERMUTATIONS P and P**T:
        //        P * inv(U**T) * inv(D) * inv(U) * P**T.
        //        Interchange rows and columns I and IPIV(I) in reverse order
        //        from the formation order of IPIV vector for Upper case.
        //
        //        ( We can use a loop over IPIV with increment 1,
        //        since the ABS value of IPIV(I) represents the row (column)
        //        index of the interchange with row (column) i in both 1x1
        //        and 2x2 pivot cases, i.e. we don't need separate code branches
        //        for 1x1 and 2x2 pivot cases )
        //
        for (i = 1; i <= n; i = i + 1) {
            ip = abs(ipiv[i - 1]);
            if (ip != i) {
                if (i < ip) {
                    Rsyswapr(uplo, n, a, lda, i, ip);
                }
                if (i > ip) {
                    Rsyswapr(uplo, n, a, lda, ip, i);
                }
            }
        }
        //
    } else {
        //
        //        Begin Lower
        //
        //        inv A = P * inv(L**T) * inv(D) * inv(L) * P**T.
        //
        Rtrtri(uplo, "U", n, a, lda, info);
        //
        //        inv(D) and inv(D) * inv(L)
        //
        k = n;
        while (k >= 1) {
            if (ipiv[k - 1] > 0) {
                //              1 x 1 diagonal NNB
                work[(k - 1) + (invd - 1) * ldwork] = one / a[(k - 1) + (k - 1) * lda];
                work[(k - 1) + ((invd + 1) - 1) * ldwork] = zero;
            } else {
                //              2 x 2 diagonal NNB
                t = work[((k - 1) - 1)];
                ak = a[((k - 1) - 1) + ((k - 1) - 1) * lda] / t;
                akp1 = a[(k - 1) + (k - 1) * lda] / t;
                akkp1 = work[((k - 1) - 1)] / t;
                d = t * (ak * akp1 - one);
                work[((k - 1) - 1) + (invd - 1) * ldwork] = akp1 / d;
                work[(k - 1) + (invd - 1) * ldwork] = ak / d;
                work[(k - 1) + ((invd + 1) - 1) * ldwork] = -akkp1 / d;
                work[((k - 1) - 1) + ((invd + 1) - 1) * ldwork] = work[(k - 1) + ((invd + 1) - 1) * ldwork];
                k = k - 1;
            }
            k = k - 1;
        }
        //
        //        inv(L**T) = (inv(L))**T
        //
        //        inv(L**T) * inv(D) * inv(L)
        //
        cut = 0;
        while (cut < n) {
            nnb = nb;
            if ((cut + nnb) > n) {
                nnb = n - cut;
            } else {
                icount = 0;
                //              count negative elements,
                for (i = cut + 1; i <= cut + nnb; i = i + 1) {
                    if (ipiv[i - 1] < 0) {
                        icount++;
                    }
                }
                //              need a even number for a clear cut
                if (mod(icount, 2) == 1) {
                    nnb++;
                }
            }
            //
            //           L21 Block
            //
            for (i = 1; i <= n - cut - nnb; i = i + 1) {
                for (j = 1; j <= nnb; j = j + 1) {
                    work[(i - 1) + (j - 1) * ldwork] = a[((cut + nnb + i) - 1) + ((cut + j) - 1) * lda];
                }
            }
            //
            //           L11 Block
            //
            for (i = 1; i <= nnb; i = i + 1) {
                work[((u11 + i) - 1) + (i - 1) * ldwork] = one;
                for (j = i + 1; j <= nnb; j = j + 1) {
                    work[((u11 + i) - 1) + (j - 1) * ldwork] = zero;
                }
                for (j = 1; j <= i - 1; j = j + 1) {
                    work[((u11 + i) - 1) + (j - 1) * ldwork] = a[((cut + i) - 1) + ((cut + j) - 1) * lda];
                }
            }
            //
            //           invD*L21
            //
            i = n - cut - nnb;
            while (i >= 1) {
                if (ipiv[(cut + nnb + i) - 1] > 0) {
                    for (j = 1; j <= nnb; j = j + 1) {
                        work[(i - 1) + (j - 1) * ldwork] = work[((cut + nnb + i) - 1) + (invd - 1) * ldwork] * work[(i - 1) + (j - 1) * ldwork];
                    }
                } else {
                    for (j = 1; j <= nnb; j = j + 1) {
                        u01_i_j = work[(i - 1) + (j - 1) * ldwork];
                        u01_ip1_j = work[((i - 1) - 1) + (j - 1) * ldwork];
                        work[(i - 1) + (j - 1) * ldwork] = work[((cut + nnb + i) - 1) + (invd - 1) * ldwork] * u01_i_j + work[((cut + nnb + i) - 1) + ((invd + 1) - 1) * ldwork] * u01_ip1_j;
                        work[((i - 1) - 1) + (j - 1) * ldwork] = work[((cut + nnb + i - 1) - 1) + ((invd + 1) - 1) * ldwork] * u01_i_j + work[((cut + nnb + i - 1) - 1) + (invd - 1) * ldwork] * u01_ip1_j;
                    }
                    i = i - 1;
                }
                i = i - 1;
            }
            //
            //           invD1*L11
            //
            i = nnb;
            while (i >= 1) {
                if (ipiv[(cut + i) - 1] > 0) {
                    for (j = 1; j <= nnb; j = j + 1) {
                        work[((u11 + i) - 1) + (j - 1) * ldwork] = work[((cut + i) - 1) + (invd - 1) * ldwork] * work[((u11 + i) - 1) + (j - 1) * ldwork];
                    }
                    //
                } else {
                    for (j = 1; j <= nnb; j = j + 1) {
                        u11_i_j = work[((u11 + i) - 1) + (j - 1) * ldwork];
                        u11_ip1_j = work[((u11 + i - 1) - 1) + (j - 1) * ldwork];
                        work[((u11 + i) - 1) + (j - 1) * ldwork] = work[((cut + i) - 1) + (invd - 1) * ldwork] * work[((u11 + i) - 1) + (j - 1) * ldwork] + work[((cut + i) - 1) + ((invd + 1) - 1) * ldwork] * u11_ip1_j;
                        work[((u11 + i - 1) - 1) + (j - 1) * ldwork] = work[((cut + i - 1) - 1) + ((invd + 1) - 1) * ldwork] * u11_i_j + work[((cut + i - 1) - 1) + (invd - 1) * ldwork] * u11_ip1_j;
                    }
                    i = i - 1;
                }
                i = i - 1;
            }
            //
            //           L11**T * invD1 * L11 -> L11
            //
            Rtrmm("L", uplo, "T", "U", nnb, nnb, one, &a[((cut + 1) - 1) + ((cut + 1) - 1) * lda], lda, &work[((u11 + 1) - 1)], n + nb + 1);
            //
            for (i = 1; i <= nnb; i = i + 1) {
                for (j = 1; j <= i; j = j + 1) {
                    a[((cut + i) - 1) + ((cut + j) - 1) * lda] = work[((u11 + i) - 1) + (j - 1) * ldwork];
                }
            }
            //
            if ((cut + nnb) < n) {
                //
                //              L21**T * invD2*L21 -> A( CUT+I, CUT+J )
                //
                Rgemm("T", "N", nnb, nnb, n - nnb - cut, one, &a[((cut + nnb + 1) - 1) + ((cut + 1) - 1) * lda], lda, work, n + nb + 1, zero, &work[((u11 + 1) - 1)], n + nb + 1);
                //
                //              L11 =  L11**T * invD1 * L11 + U01**T * invD * U01
                //
                for (i = 1; i <= nnb; i = i + 1) {
                    for (j = 1; j <= i; j = j + 1) {
                        a[((cut + i) - 1) + ((cut + j) - 1) * lda] += work[((u11 + i) - 1) + (j - 1) * ldwork];
                    }
                }
                //
                //              L01 =  L22**T * invD2 * L21
                //
                Rtrmm("L", uplo, "T", "U", n - nnb - cut, nnb, one, &a[((cut + nnb + 1) - 1) + ((cut + nnb + 1) - 1) * lda], lda, work, n + nb + 1);
                //
                //              Update L21
                //
                for (i = 1; i <= n - cut - nnb; i = i + 1) {
                    for (j = 1; j <= nnb; j = j + 1) {
                        a[((cut + nnb + i) - 1) + ((cut + j) - 1) * lda] = work[(i - 1) + (j - 1) * ldwork];
                    }
                }
                //
            } else {
                //
                //              L11 =  L11**T * invD1 * L11
                //
                for (i = 1; i <= nnb; i = i + 1) {
                    for (j = 1; j <= i; j = j + 1) {
                        a[((cut + i) - 1) + ((cut + j) - 1) * lda] = work[((u11 + i) - 1) + (j - 1) * ldwork];
                    }
                }
            }
            //
            //           Next Block
            //
            cut += nnb;
            //
        }
        //
        //        Apply PERMUTATIONS P and P**T:
        //        P * inv(L**T) * inv(D) * inv(L) * P**T.
        //        Interchange rows and columns I and IPIV(I) in reverse order
        //        from the formation order of IPIV vector for Lower case.
        //
        //        ( We can use a loop over IPIV with increment -1,
        //        since the ABS value of IPIV(I) represents the row (column)
        //        index of the interchange with row (column) i in both 1x1
        //        and 2x2 pivot cases, i.e. we don't need separate code branches
        //        for 1x1 and 2x2 pivot cases )
        //
        for (i = n; i >= 1; i = i - 1) {
            ip = abs(ipiv[i - 1]);
            if (ip != i) {
                if (i < ip) {
                    Rsyswapr(uplo, n, a, lda, i, ip);
                }
                if (i > ip) {
                    Rsyswapr(uplo, n, a, lda, ip, i);
                }
            }
        }
        //
    }
    //
    //     End of Rsytri_3x
    //
}

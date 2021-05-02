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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Rlavsp(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER *ipiv, REAL *b, INTEGER const ldb, INTEGER &info) {
    bool nounit = false;
    INTEGER k = 0;
    INTEGER kc = 0;
    const REAL one = 1.0;
    INTEGER kp = 0;
    INTEGER kcnext = 0;
    REAL d11 = 0.0;
    REAL d22 = 0.0;
    REAL d12 = 0.0;
    REAL d21 = 0.0;
    INTEGER j = 0;
    REAL t1 = 0.0;
    REAL t2 = 0.0;
    //
    //  -- LAPACK test routine --
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
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (!Mlsame(trans, "N") && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = -2;
    } else if (!Mlsame(diag, "U") && !Mlsame(diag, "N")) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Rlavsp ", -info);
        return;
    }
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        return;
    }
    //
    nounit = Mlsame(diag, "N");
    //------------------------------------------
    //
    //     Compute  B := A * B  (No transpose)
    //
    //------------------------------------------
    if (Mlsame(trans, "N")) {
        //
        //        Compute  B := U*B
        //        where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
        //
        if (Mlsame(uplo, "U")) {
            //
            //        Loop forward applying the transformations.
            //
            k = 1;
            kc = 1;
        statement_10:
            if (k > n) {
                goto statement_30;
            }
            //
            //           1 x 1 pivot block
            //
            if (ipiv[k - 1] > 0) {
                //
                //              Multiply by the diagonal element if forming U * D.
                //
                if (nounit) {
                    Rscal(nrhs, &a[(kc + k - 1) - 1], &b[(k - 1)], ldb);
                }
                //
                //              Multiply by P(K) * inv(U(K))  if K > 1.
                //
                if (k > 1) {
                    //
                    //                 Apply the transformation.
                    //
                    Rger(k - 1, nrhs, one, &a[kc - 1], 1, &b[(k - 1)], ldb, &b[(1 - 1)], ldb);
                    //
                    //                 Interchange if P(K) != I.
                    //
                    kp = ipiv[k - 1];
                    if (kp != k) {
                        Rswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                }
                kc += k;
                k++;
            } else {
                //
                //              2 x 2 pivot block
                //
                kcnext = kc + k;
                //
                //              Multiply by the diagonal block if forming U * D.
                //
                if (nounit) {
                    d11 = a[(kcnext - 1) - 1];
                    d22 = a[(kcnext + k) - 1];
                    d12 = a[(kcnext + k - 1) - 1];
                    d21 = d12;
                    for (j = 1; j <= nrhs; j = j + 1) {
                        t1 = b[(k - 1) + (j - 1) * ldb];
                        t2 = b[((k + 1) - 1) + (j - 1) * ldb];
                        b[(k - 1) + (j - 1) * ldb] = d11 * t1 + d12 * t2;
                        b[((k + 1) - 1) + (j - 1) * ldb] = d21 * t1 + d22 * t2;
                    }
                }
                //
                //              Multiply by  P(K) * inv(U(K))  if K > 1.
                //
                if (k > 1) {
                    //
                    //                 Apply the transformations.
                    //
                    Rger(k - 1, nrhs, one, &a[kc - 1], 1, &b[(k - 1)], ldb, &b[(1 - 1)], ldb);
                    Rger(k - 1, nrhs, one, &a[kcnext - 1], 1, &b[((k + 1) - 1)], ldb, &b[(1 - 1)], ldb);
                    //
                    //                 Interchange if P(K) != I.
                    //
                    kp = abs(ipiv[k - 1]);
                    if (kp != k) {
                        Rswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                }
                kc = kcnext + k + 1;
                k += 2;
            }
            goto statement_10;
        statement_30:;
            //
            //        Compute  B := L*B
            //        where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
            //
        } else {
            //
            //           Loop backward applying the transformations to B.
            //
            k = n;
            kc = n * (n + 1) / 2 + 1;
        statement_40:
            if (k < 1) {
                goto statement_60;
            }
            kc = kc - (n - k + 1);
            //
            //           Test the pivot index.  If greater than zero, a 1 x 1
            //           pivot was used, otherwise a 2 x 2 pivot was used.
            //
            if (ipiv[k - 1] > 0) {
                //
                //              1 x 1 pivot block:
                //
                //              Multiply by the diagonal element if forming L * D.
                //
                if (nounit) {
                    Rscal(nrhs, &a[kc - 1], &b[(k - 1)], ldb);
                }
                //
                //              Multiply by  P(K) * inv(L(K))  if K < N.
                //
                if (k != n) {
                    kp = ipiv[k - 1];
                    //
                    //                 Apply the transformation.
                    //
                    Rger(n - k, nrhs, one, &a[(kc + 1) - 1], 1, &b[(k - 1)], ldb, &b[((k + 1) - 1)], ldb);
                    //
                    //                 Interchange if a permutation was applied at the
                    //                 K-th step of the factorization.
                    //
                    if (kp != k) {
                        Rswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                }
                k = k - 1;
                //
            } else {
                //
                //              2 x 2 pivot block:
                //
                kcnext = kc - (n - k + 2);
                //
                //              Multiply by the diagonal block if forming L * D.
                //
                if (nounit) {
                    d11 = a[kcnext - 1];
                    d22 = a[kc - 1];
                    d21 = a[(kcnext + 1) - 1];
                    d12 = d21;
                    for (j = 1; j <= nrhs; j = j + 1) {
                        t1 = b[((k - 1) - 1) + (j - 1) * ldb];
                        t2 = b[(k - 1) + (j - 1) * ldb];
                        b[((k - 1) - 1) + (j - 1) * ldb] = d11 * t1 + d12 * t2;
                        b[(k - 1) + (j - 1) * ldb] = d21 * t1 + d22 * t2;
                    }
                }
                //
                //              Multiply by  P(K) * inv(L(K))  if K < N.
                //
                if (k != n) {
                    //
                    //                 Apply the transformation.
                    //
                    Rger(n - k, nrhs, one, &a[(kc + 1) - 1], 1, &b[(k - 1)], ldb, &b[((k + 1) - 1)], ldb);
                    Rger(n - k, nrhs, one, &a[(kcnext + 2) - 1], 1, &b[((k - 1) - 1)], ldb, &b[((k + 1) - 1)], ldb);
                    //
                    //                 Interchange if a permutation was applied at the
                    //                 K-th step of the factorization.
                    //
                    kp = abs(ipiv[k - 1]);
                    if (kp != k) {
                        Rswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                }
                kc = kcnext;
                k = k - 2;
            }
            goto statement_40;
        statement_60:;
        }
        //----------------------------------------
        //
        //     Compute  B := A' * B  (transpose)
        //
        //----------------------------------------
    } else {
        //
        //        Form  B := U'*B
        //        where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
        //        and   U' = inv(U'(1))*P(1)* ... *inv(U'(m))*P(m)
        //
        if (Mlsame(uplo, "U")) {
            //
            //           Loop backward applying the transformations.
            //
            k = n;
            kc = n * (n + 1) / 2 + 1;
        statement_70:
            if (k < 1) {
                goto statement_90;
            }
            kc = kc - k;
            //
            //           1 x 1 pivot block.
            //
            if (ipiv[k - 1] > 0) {
                if (k > 1) {
                    //
                    //                 Interchange if P(K) != I.
                    //
                    kp = ipiv[k - 1];
                    if (kp != k) {
                        Rswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                    //
                    //                 Apply the transformation
                    //
                    Rgemv("Transpose", k - 1, nrhs, one, b, ldb, &a[kc - 1], 1, one, &b[(k - 1)], ldb);
                }
                if (nounit) {
                    Rscal(nrhs, &a[(kc + k - 1) - 1], &b[(k - 1)], ldb);
                }
                k = k - 1;
                //
                //           2 x 2 pivot block.
                //
            } else {
                kcnext = kc - (k - 1);
                if (k > 2) {
                    //
                    //                 Interchange if P(K) != I.
                    //
                    kp = abs(ipiv[k - 1]);
                    if (kp != k - 1) {
                        Rswap(nrhs, &b[((k - 1) - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                    //
                    //                 Apply the transformations
                    //
                    Rgemv("Transpose", k - 2, nrhs, one, b, ldb, &a[kc - 1], 1, one, &b[(k - 1)], ldb);
                    Rgemv("Transpose", k - 2, nrhs, one, b, ldb, &a[kcnext - 1], 1, one, &b[((k - 1) - 1)], ldb);
                }
                //
                //              Multiply by the diagonal block if non-unit.
                //
                if (nounit) {
                    d11 = a[(kc - 1) - 1];
                    d22 = a[(kc + k - 1) - 1];
                    d12 = a[(kc + k - 2) - 1];
                    d21 = d12;
                    for (j = 1; j <= nrhs; j = j + 1) {
                        t1 = b[((k - 1) - 1) + (j - 1) * ldb];
                        t2 = b[(k - 1) + (j - 1) * ldb];
                        b[((k - 1) - 1) + (j - 1) * ldb] = d11 * t1 + d12 * t2;
                        b[(k - 1) + (j - 1) * ldb] = d21 * t1 + d22 * t2;
                    }
                }
                kc = kcnext;
                k = k - 2;
            }
            goto statement_70;
        statement_90:;
            //
            //        Form  B := L'*B
            //        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
            //        and   L' = inv(L(m))*P(m)* ... *inv(L(1))*P(1)
            //
        } else {
            //
            //           Loop forward applying the L-transformations.
            //
            k = 1;
            kc = 1;
        statement_100:
            if (k > n) {
                goto statement_120;
            }
            //
            //           1 x 1 pivot block
            //
            if (ipiv[k - 1] > 0) {
                if (k < n) {
                    //
                    //                 Interchange if P(K) != I.
                    //
                    kp = ipiv[k - 1];
                    if (kp != k) {
                        Rswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                    //
                    //                 Apply the transformation
                    //
                    Rgemv("Transpose", n - k, nrhs, one, &b[((k + 1) - 1)], ldb, &a[(kc + 1) - 1], 1, one, &b[(k - 1)], ldb);
                }
                if (nounit) {
                    Rscal(nrhs, &a[kc - 1], &b[(k - 1)], ldb);
                }
                kc += n - k + 1;
                k++;
                //
                //           2 x 2 pivot block.
                //
            } else {
                kcnext = kc + n - k + 1;
                if (k < n - 1) {
                    //
                    //              Interchange if P(K) != I.
                    //
                    kp = abs(ipiv[k - 1]);
                    if (kp != k + 1) {
                        Rswap(nrhs, &b[((k + 1) - 1)], ldb, &b[(kp - 1)], ldb);
                    }
                    //
                    //                 Apply the transformation
                    //
                    Rgemv("Transpose", n - k - 1, nrhs, one, &b[((k + 2) - 1)], ldb, &a[(kcnext + 1) - 1], 1, one, &b[((k + 1) - 1)], ldb);
                    Rgemv("Transpose", n - k - 1, nrhs, one, &b[((k + 2) - 1)], ldb, &a[(kc + 2) - 1], 1, one, &b[(k - 1)], ldb);
                }
                //
                //              Multiply by the diagonal block if non-unit.
                //
                if (nounit) {
                    d11 = a[kc - 1];
                    d22 = a[kcnext - 1];
                    d21 = a[(kc + 1) - 1];
                    d12 = d21;
                    for (j = 1; j <= nrhs; j = j + 1) {
                        t1 = b[(k - 1) + (j - 1) * ldb];
                        t2 = b[((k + 1) - 1) + (j - 1) * ldb];
                        b[(k - 1) + (j - 1) * ldb] = d11 * t1 + d12 * t2;
                        b[((k + 1) - 1) + (j - 1) * ldb] = d21 * t1 + d22 * t2;
                    }
                }
                kc = kcnext + (n - k);
                k += 2;
            }
            goto statement_100;
        statement_120:;
        }
        //
    }
    //
    //     End of Rlavsp
    //
}

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

void Chpgst(INTEGER const itype, const char *uplo, INTEGER const n, COMPLEX *ap, COMPLEX *bp, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    if (itype < 1 || itype > 3) {
        info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    }
    if (info != 0) {
        Mxerbla("Chpgst", -info);
        return;
    }
    //
    INTEGER jj = 0;
    INTEGER j = 0;
    INTEGER j1 = 0;
    REAL bjj = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const REAL one = 1.0;
    INTEGER kk = 0;
    INTEGER k = 0;
    INTEGER k1k1 = 0;
    REAL akk = 0.0;
    REAL bkk = 0.0;
    const REAL half = 0.5e+0;
    COMPLEX ct = 0.0;
    INTEGER k1 = 0;
    INTEGER j1j1 = 0;
    REAL ajj = 0.0;
    if (itype == 1) {
        if (upper) {
            //
            //           Compute inv(U**H)*A*inv(U)
            //
            //           J1 and JJ are the indices of A(1,j) and A(j,j)
            //
            jj = 0;
            for (j = 1; j <= n; j = j + 1) {
                j1 = jj + 1;
                jj += j;
                //
                //              Compute the j-th column of the upper triangle of A
                //
                ap[jj - 1] = ap[jj - 1].real();
                bjj = bp[jj - 1].real();
                Ctpsv(uplo, "Conjugate transpose", "Non-unit", j, bp, &ap[j1 - 1], 1);
                Chpmv(uplo, j - 1, -cone, ap, &bp[j1 - 1], 1, cone, &ap[j1 - 1], 1);
                CRscal(j - 1, one / bjj, &ap[j1 - 1], 1);
                ap[jj - 1] = (ap[jj - 1] - Cdotc(j - 1, &ap[j1 - 1], 1, &bp[j1 - 1], 1)) / bjj;
            }
        } else {
            //
            //           Compute inv(L)*A*inv(L**H)
            //
            //           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)
            //
            kk = 1;
            for (k = 1; k <= n; k = k + 1) {
                k1k1 = kk + n - k + 1;
                //
                //              Update the lower triangle of A(k:n,k:n)
                //
                akk = ap[kk - 1].real();
                bkk = bp[kk - 1].real();
                akk = akk / pow2(bkk);
                ap[kk - 1] = akk;
                if (k < n) {
                    CRscal(n - k, one / bkk, &ap[(kk + 1) - 1], 1);
                    ct = -half * akk;
                    Caxpy(n - k, ct, &bp[(kk + 1) - 1], 1, &ap[(kk + 1) - 1], 1);
                    Chpr2(uplo, n - k, -cone, &ap[(kk + 1) - 1], 1, &bp[(kk + 1) - 1], 1, &ap[k1k1 - 1]);
                    Caxpy(n - k, ct, &bp[(kk + 1) - 1], 1, &ap[(kk + 1) - 1], 1);
                    Ctpsv(uplo, "No transpose", "Non-unit", n - k, &bp[k1k1 - 1], &ap[(kk + 1) - 1], 1);
                }
                kk = k1k1;
            }
        }
    } else {
        if (upper) {
            //
            //           Compute U*A*U**H
            //
            //           K1 and KK are the indices of A(1,k) and A(k,k)
            //
            kk = 0;
            for (k = 1; k <= n; k = k + 1) {
                k1 = kk + 1;
                kk += k;
                //
                //              Update the upper triangle of A(1:k,1:k)
                //
                akk = ap[kk - 1].real();
                bkk = bp[kk - 1].real();
                Ctpmv(uplo, "No transpose", "Non-unit", k - 1, bp, &ap[k1 - 1], 1);
                ct = half * akk;
                Caxpy(k - 1, ct, &bp[k1 - 1], 1, &ap[k1 - 1], 1);
                Chpr2(uplo, k - 1, cone, &ap[k1 - 1], 1, &bp[k1 - 1], 1, ap);
                Caxpy(k - 1, ct, &bp[k1 - 1], 1, &ap[k1 - 1], 1);
                CRscal(k - 1, bkk, &ap[k1 - 1], 1);
                ap[kk - 1] = akk * pow2(bkk);
            }
        } else {
            //
            //           Compute L**H *A*L
            //
            //           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)
            //
            jj = 1;
            for (j = 1; j <= n; j = j + 1) {
                j1j1 = jj + n - j + 1;
                //
                //              Compute the j-th column of the lower triangle of A
                //
                ajj = ap[jj - 1].real();
                bjj = bp[jj - 1].real();
                ap[jj - 1] = ajj * bjj + Cdotc(n - j, &ap[(jj + 1) - 1], 1, &bp[(jj + 1) - 1], 1);
                CRscal(n - j, bjj, &ap[(jj + 1) - 1], 1);
                Chpmv(uplo, n - j, cone, &ap[j1j1 - 1], &bp[(jj + 1) - 1], 1, cone, &ap[(jj + 1) - 1], 1);
                Ctpmv(uplo, "Conjugate transpose", "Non-unit", n - j + 1, &bp[jj - 1], &ap[jj - 1], 1);
                jj = j1j1;
            }
        }
    }
    //
    //     End of Chpgst
    //
}

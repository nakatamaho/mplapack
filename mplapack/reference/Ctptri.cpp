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

void Ctptri(const char *uplo, const char *diag, INTEGER const &n, COMPLEX *ap, INTEGER &info) {
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
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    bool nounit = Mlsame(diag, "N");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (!nounit && !Mlsame(diag, "U")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    }
    if (info != 0) {
        Mxerbla("Ctptri", -info);
        return;
    }
    //
    //     Check for singularity if non-unit.
    //
    INTEGER jj = 0;
    const COMPLEX zero = (0.0, 0.0);
    if (nounit) {
        if (upper) {
            jj = 0;
            for (info = 1; info <= n; info = info + 1) {
                jj += info;
                if (ap[jj - 1] == zero) {
                    return;
                }
            }
        } else {
            jj = 1;
            for (info = 1; info <= n; info = info + 1) {
                if (ap[jj - 1] == zero) {
                    return;
                }
                jj += n - info + 1;
            }
        }
        info = 0;
    }
    //
    INTEGER jc = 0;
    INTEGER j = 0;
    const COMPLEX one = (1.0, 0.0);
    COMPLEX ajj = 0.0;
    INTEGER jclast = 0;
    if (upper) {
        //
        //        Compute inverse of upper triangular matrix.
        //
        jc = 1;
        for (j = 1; j <= n; j = j + 1) {
            if (nounit) {
                ap[(jc + j - 1) - 1] = one / ap[(jc + j - 1) - 1];
                ajj = -ap[(jc + j - 1) - 1];
            } else {
                ajj = -one;
            }
            //
            //           Compute elements 1:j-1 of j-th column.
            //
            Ctpmv("Upper", "No transpose", diag, j - 1, ap, ap[jc - 1], 1);
            Cscal(j - 1, ajj, ap[jc - 1], 1);
            jc += j;
        }
        //
    } else {
        //
        //        Compute inverse of lower triangular matrix.
        //
        jc = n * (n + 1) / 2;
        for (j = n; j >= 1; j = j - 1) {
            if (nounit) {
                ap[jc - 1] = one / ap[jc - 1];
                ajj = -ap[jc - 1];
            } else {
                ajj = -one;
            }
            if (j < n) {
                //
                //              Compute elements j+1:n of j-th column.
                //
                Ctpmv("Lower", "No transpose", diag, n - j, ap[jclast - 1], ap[(jc + 1) - 1], 1);
                Cscal(n - j, ajj, ap[(jc + 1) - 1], 1);
            }
            jclast = jc;
            jc = jc - n + j - 2;
        }
    }
    //
    //     End of Ctptri
    //
}

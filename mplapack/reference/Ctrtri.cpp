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

void Ctrtri(const char *uplo, const char *diag, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER &info) {
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
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Ctrtri", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Check for singularity if non-unit.
    //
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    if (nounit) {
        for (info = 1; info <= n; info = info + 1) {
            if (a[(info - 1) + (info - 1) * lda] == zero) {
                return;
            }
        }
        info = 0;
    }
    //
    //     Determine the block size for this environment.
    //
    char uplo_diag[3];
    uplo_diag[0] = uplo[0];
    uplo_diag[1] = diag[0];
    uplo_diag[2] = '\0';
    INTEGER nb = iMlaenv(1, "Ctrtri", uplo_diag, n, -1, -1, -1);
    INTEGER j = 0;
    INTEGER jb = 0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    INTEGER nn = 0;
    if (nb <= 1 || nb >= n) {
        //
        //        Use unblocked code
        //
        Ctrti2(uplo, diag, n, a, lda, info);
    } else {
        //
        //        Use blocked code
        //
        if (upper) {
            //
            //           Compute inverse of upper triangular matrix
            //
            for (j = 1; j <= n; j = j + nb) {
                jb = min(nb, n - j + 1);
                //
                //              Compute rows 1:j-1 of current block column
                //
                Ctrmm("Left", "Upper", "No transpose", diag, j - 1, jb, one, a, lda, &a[(j - 1) * lda], lda);
                Ctrsm("Right", "Upper", "No transpose", diag, j - 1, jb, -one, &a[(j - 1) + (j - 1) * lda], lda, &a[(j - 1) * lda], lda);
                //
                //              Compute inverse of current diagonal block
                //
                Ctrti2("Upper", diag, jb, &a[(j - 1) + (j - 1) * lda], lda, info);
            }
        } else {
            //
            //           Compute inverse of lower triangular matrix
            //
            nn = ((n - 1) / nb) * nb + 1;
            for (j = nn; j >= 1; j = j - nb) {
                jb = min(nb, n - j + 1);
                if (j + jb <= n) {
                    //
                    //                 Compute rows j+jb:n of current block column
                    //
                    Ctrmm("Left", "Lower", "No transpose", diag, n - j - jb + 1, jb, one, &a[((j + jb) - 1) + ((j + jb) - 1) * lda], lda, &a[((j + jb) - 1) + (j - 1) * lda], lda);
                    Ctrsm("Right", "Lower", "No transpose", diag, n - j - jb + 1, jb, -one, &a[(j - 1) + (j - 1) * lda], lda, &a[((j + jb) - 1) + (j - 1) * lda], lda);
                }
                //
                //              Compute inverse of current diagonal block
                //
                Ctrti2("Lower", diag, jb, &a[(j - 1) + (j - 1) * lda], lda, info);
            }
        }
    }
    //
    //     End of Ctrtri
    //
}

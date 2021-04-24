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

void Clakf2(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *d, COMPLEX *e, COMPLEX *z, INTEGER const ldz) {
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
    //  ====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize Z
    //
    INTEGER mn = m * n;
    INTEGER mn2 = 2 * mn;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    zlaset("Full", mn2, mn2, zero, zero, z, ldz);
    //
    INTEGER ik = 1;
    INTEGER l = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    for (l = 1; l <= n; l = l + 1) {
        //
        //        form kron(In, A)
        //
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= m; j = j + 1) {
                z[((ik + i - 1) - 1) + ((ik + j - 1) - 1) * ldz] = a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        //        form kron(In, D)
        //
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= m; j = j + 1) {
                z[((ik + mn + i - 1) - 1) + ((ik + j - 1) - 1) * ldz] = d[(i - 1) + (j - 1) * ldd];
            }
        }
        //
        ik += m;
    }
    //
    ik = 1;
    INTEGER jk = 0;
    for (l = 1; l <= n; l = l + 1) {
        jk = mn + 1;
        //
        for (j = 1; j <= n; j = j + 1) {
            //
            //           form -kron(B', Im)
            //
            for (i = 1; i <= m; i = i + 1) {
                z[((ik + i - 1) - 1) + ((jk + i - 1) - 1) * ldz] = -b[(j - 1) + (l - 1) * ldb];
            }
            //
            //           form -kron(E', Im)
            //
            for (i = 1; i <= m; i = i + 1) {
                z[((ik + mn + i - 1) - 1) + ((jk + i - 1) - 1) * ldz] = -e[(j - 1) + (l - 1) * lde];
            }
            //
            jk += m;
        }
        //
        ik += m;
    }
    //
    //     End of Clakf2
    //
}

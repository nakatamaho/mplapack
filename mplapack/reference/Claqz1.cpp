/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZLAQZ1.
// Original LAPACK authors:
//   Thijs Steel, KU Leuven

#include <mpblas.h>
#include <mplapack.h>

void Claqz1(bool const ilq, bool const ilz, INTEGER const k, INTEGER const istartm, INTEGER const istopm, INTEGER const ihi, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, INTEGER const nq, INTEGER const qstart, COMPLEX *q, INTEGER const ldq, INTEGER const nz, INTEGER const zstart, COMPLEX *z, INTEGER const ldz) {
    //
    // Arguments
    //
    // Parameters
    //
    // Local variables
    //
    // External Functions
    //
    REAL c = 0.0;
    COMPLEX s = 0.0;
    COMPLEX temp = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (k + 1 == ihi) {
        //
        // Shift is located on the edge of the matrix, remove it
        //
        Clartg(b[(ihi - 1) + (ihi - 1) * ldb], b[(ihi - 1) + ((ihi - 1) - 1) * ldb], c, s, temp);
        b[(ihi - 1) + (ihi - 1) * ldb] = temp;
        b[(ihi - 1) + ((ihi - 1) - 1) * ldb] = czero;
        Crot(ihi - istartm, &b[(istartm - 1) + (ihi - 1) * ldb], 1, &b[(istartm - 1) + ((ihi - 1) - 1) * ldb], 1, c, s);
        Crot(ihi - istartm + 1, &a[(istartm - 1) + (ihi - 1) * lda], 1, &a[(istartm - 1) + ((ihi - 1) - 1) * lda], 1, c, s);
        if (ilz) {
            Crot(nz, &z[((ihi - zstart + 1) - 1) * ldz], 1, &z[((ihi - 1 - zstart + 1) - 1) * ldz], 1, c, s);
        }
        //
    } else {
        //
        // Normal operation, move bulge down
        //
        // Apply transformation from the right
        //
        Clartg(b[((k + 1) - 1) + ((k + 1) - 1) * ldb], b[((k + 1) - 1) + (k - 1) * ldb], c, s, temp);
        b[((k + 1) - 1) + ((k + 1) - 1) * ldb] = temp;
        b[((k + 1) - 1) + (k - 1) * ldb] = czero;
        Crot(k + 2 - istartm + 1, &a[(istartm - 1) + ((k + 1) - 1) * lda], 1, &a[(istartm - 1) + (k - 1) * lda], 1, c, s);
        Crot(k - istartm + 1, &b[(istartm - 1) + ((k + 1) - 1) * ldb], 1, &b[(istartm - 1) + (k - 1) * ldb], 1, c, s);
        if (ilz) {
            Crot(nz, &z[((k + 1 - zstart + 1) - 1) * ldz], 1, &z[((k - zstart + 1) - 1) * ldz], 1, c, s);
        }
        //
        // Apply transformation from the left
        //
        Clartg(a[((k + 1) - 1) + (k - 1) * lda], a[((k + 2) - 1) + (k - 1) * lda], c, s, temp);
        a[((k + 1) - 1) + (k - 1) * lda] = temp;
        a[((k + 2) - 1) + (k - 1) * lda] = czero;
        Crot(istopm - k, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda, &a[((k + 2) - 1) + ((k + 1) - 1) * lda], lda, c, s);
        Crot(istopm - k, &b[((k + 1) - 1) + ((k + 1) - 1) * ldb], ldb, &b[((k + 2) - 1) + ((k + 1) - 1) * ldb], ldb, c, s);
        if (ilq) {
            Crot(nq, &q[((k + 1 - qstart + 1) - 1) * ldq], 1, &q[((k + 2 - qstart + 1) - 1) * ldq], 1, c, conj(s));
        }
        //
    }
    //
    // End of Claqz1
    //
}

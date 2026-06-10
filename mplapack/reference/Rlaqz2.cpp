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

// Derived from LAPACK routine DLAQZ2.
// Original LAPACK authors:
//   Thijs Steel, KU Leuven

#include <mpblas.h>
#include <mplapack.h>

void Rlaqz2(bool const ilq, bool const ilz, INTEGER const k, INTEGER const istartm, INTEGER const istopm, INTEGER const ihi, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, INTEGER const nq, INTEGER const qstart, REAL *q, INTEGER const ldq, INTEGER const nz, INTEGER const zstart, REAL *z, INTEGER const ldz) {
    INTEGER ldh = 2;
    REAL h[2 * 3];
    REAL c1 = 0.0;
    REAL s1 = 0.0;
    REAL temp = 0.0;
    const REAL zero = 0.0;
    REAL c2 = 0.0;
    REAL s2 = 0.0;
    if (k + 2 == ihi) {
        // Shift is located on the edge of the matrix, remove it
        // H = B( IHI-1:IHI, IHI-2:IHI )  -- 2x3 submatrix
        h[(1 - 1) + (1 - 1) * 2] = b[((ihi - 1) - 1) + ((ihi - 2) - 1) * ldb];
        h[(2 - 1) + (1 - 1) * 2] = b[((ihi)-1) + ((ihi - 2) - 1) * ldb];
        h[(1 - 1) + (2 - 1) * 2] = b[((ihi - 1) - 1) + ((ihi - 1) - 1) * ldb];
        h[(2 - 1) + (2 - 1) * 2] = b[((ihi)-1) + ((ihi - 1) - 1) * ldb];
        h[(1 - 1) + (3 - 1) * 2] = b[((ihi - 1) - 1) + ((ihi)-1) * ldb];
        h[(2 - 1) + (3 - 1) * 2] = b[((ihi)-1) + ((ihi)-1) * ldb];
        // Make H upper triangular
        Rlartg(h[0], h[(2 - 1)], c1, s1, temp);
        h[(2 - 1)] = zero;
        h[0] = temp;
        Rrot(2, &h[(2 - 1) * ldh], 2, &h[(2 - 1) + (2 - 1) * ldh], 2, c1, s1);
        //
        Rlartg(h[(2 - 1) + (3 - 1) * ldh], h[(2 - 1) + (2 - 1) * ldh], c1, s1, temp);
        Rrot(1, &h[(3 - 1) * ldh], 1, &h[(2 - 1) * ldh], 1, c1, s1);
        Rlartg(h[(2 - 1) * ldh], h[0], c2, s2, temp);
        //
        Rrot(ihi - istartm + 1, &b[(istartm - 1) + (ihi - 1) * ldb], 1, &b[(istartm - 1) + ((ihi - 1) - 1) * ldb], 1, c1, s1);
        Rrot(ihi - istartm + 1, &b[(istartm - 1) + ((ihi - 1) - 1) * ldb], 1, &b[(istartm - 1) + ((ihi - 2) - 1) * ldb], 1, c2, s2);
        b[((ihi - 1) - 1) + ((ihi - 2) - 1) * ldb] = zero;
        b[(ihi - 1) + ((ihi - 2) - 1) * ldb] = zero;
        Rrot(ihi - istartm + 1, &a[(istartm - 1) + (ihi - 1) * lda], 1, &a[(istartm - 1) + ((ihi - 1) - 1) * lda], 1, c1, s1);
        Rrot(ihi - istartm + 1, &a[(istartm - 1) + ((ihi - 1) - 1) * lda], 1, &a[(istartm - 1) + ((ihi - 2) - 1) * lda], 1, c2, s2);
        if (ilz) {
            Rrot(nz, &z[((ihi - zstart + 1) - 1) * ldz], 1, &z[((ihi - 1 - zstart + 1) - 1) * ldz], 1, c1, s1);
            Rrot(nz, &z[((ihi - 1 - zstart + 1) - 1) * ldz], 1, &z[((ihi - 2 - zstart + 1) - 1) * ldz], 1, c2, s2);
        }
        //
        Rlartg(a[((ihi - 1) - 1) + ((ihi - 2) - 1) * lda], a[(ihi - 1) + ((ihi - 2) - 1) * lda], c1, s1, temp);
        a[((ihi - 1) - 1) + ((ihi - 2) - 1) * lda] = temp;
        a[(ihi - 1) + ((ihi - 2) - 1) * lda] = zero;
        Rrot(istopm - ihi + 2, &a[((ihi - 1) - 1) + ((ihi - 1) - 1) * lda], lda, &a[(ihi - 1) + ((ihi - 1) - 1) * lda], lda, c1, s1);
        Rrot(istopm - ihi + 2, &b[((ihi - 1) - 1) + ((ihi - 1) - 1) * ldb], ldb, &b[(ihi - 1) + ((ihi - 1) - 1) * ldb], ldb, c1, s1);
        if (ilq) {
            Rrot(nq, &q[((ihi - 1 - qstart + 1) - 1) * ldq], 1, &q[((ihi - qstart + 1) - 1) * ldq], 1, c1, s1);
        }
        //
        Rlartg(b[(ihi - 1) + (ihi - 1) * ldb], b[(ihi - 1) + ((ihi - 1) - 1) * ldb], c1, s1, temp);
        b[(ihi - 1) + (ihi - 1) * ldb] = temp;
        b[(ihi - 1) + ((ihi - 1) - 1) * ldb] = zero;
        Rrot(ihi - istartm, &b[(istartm - 1) + (ihi - 1) * ldb], 1, &b[(istartm - 1) + ((ihi - 1) - 1) * ldb], 1, c1, s1);
        Rrot(ihi - istartm + 1, &a[(istartm - 1) + (ihi - 1) * lda], 1, &a[(istartm - 1) + ((ihi - 1) - 1) * lda], 1, c1, s1);
        if (ilz) {
            Rrot(nz, &z[((ihi - zstart + 1) - 1) * ldz], 1, &z[((ihi - 1 - zstart + 1) - 1) * ldz], 1, c1, s1);
        }
        //
    } else {
        //
        // Normal operation, move bulge down
        //
        // H = B( K+1:K+2, K:K+2 )  -- 2x3 submatrix
        h[(1 - 1) + (1 - 1) * 2] = b[((k + 1) - 1) + ((k)-1) * ldb];
        h[(2 - 1) + (1 - 1) * 2] = b[((k + 2) - 1) + ((k)-1) * ldb];
        h[(1 - 1) + (2 - 1) * 2] = b[((k + 1) - 1) + ((k + 1) - 1) * ldb];
        h[(2 - 1) + (2 - 1) * 2] = b[((k + 2) - 1) + ((k + 1) - 1) * ldb];
        h[(1 - 1) + (3 - 1) * 2] = b[((k + 1) - 1) + ((k + 2) - 1) * ldb];
        h[(2 - 1) + (3 - 1) * 2] = b[((k + 2) - 1) + ((k + 2) - 1) * ldb];
        //
        // Make H upper triangular
        //
        Rlartg(h[0], h[(2 - 1)], c1, s1, temp);
        h[(2 - 1)] = zero;
        h[0] = temp;
        Rrot(2, &h[(2 - 1) * ldh], 2, &h[(2 - 1) + (2 - 1) * ldh], 2, c1, s1);
        //
        // Calculate Z1 and Z2
        //
        Rlartg(h[(2 - 1) + (3 - 1) * ldh], h[(2 - 1) + (2 - 1) * ldh], c1, s1, temp);
        Rrot(1, &h[(3 - 1) * ldh], 1, &h[(2 - 1) * ldh], 1, c1, s1);
        Rlartg(h[(2 - 1) * ldh], h[0], c2, s2, temp);
        //
        // Apply transformations from the right
        //
        Rrot(k + 3 - istartm + 1, &a[(istartm - 1) + ((k + 2) - 1) * lda], 1, &a[(istartm - 1) + ((k + 1) - 1) * lda], 1, c1, s1);
        Rrot(k + 3 - istartm + 1, &a[(istartm - 1) + ((k + 1) - 1) * lda], 1, &a[(istartm - 1) + (k - 1) * lda], 1, c2, s2);
        Rrot(k + 2 - istartm + 1, &b[(istartm - 1) + ((k + 2) - 1) * ldb], 1, &b[(istartm - 1) + ((k + 1) - 1) * ldb], 1, c1, s1);
        Rrot(k + 2 - istartm + 1, &b[(istartm - 1) + ((k + 1) - 1) * ldb], 1, &b[(istartm - 1) + (k - 1) * ldb], 1, c2, s2);
        if (ilz) {
            Rrot(nz, &z[((k + 2 - zstart + 1) - 1) * ldz], 1, &z[((k + 1 - zstart + 1) - 1) * ldz], 1, c1, s1);
            Rrot(nz, &z[((k + 1 - zstart + 1) - 1) * ldz], 1, &z[((k - zstart + 1) - 1) * ldz], 1, c2, s2);
        }
        b[((k + 1) - 1) + (k - 1) * ldb] = zero;
        b[((k + 2) - 1) + (k - 1) * ldb] = zero;
        //
        // Calculate Q1 and Q2
        //
        Rlartg(a[((k + 2) - 1) + (k - 1) * lda], a[((k + 3) - 1) + (k - 1) * lda], c1, s1, temp);
        a[((k + 2) - 1) + (k - 1) * lda] = temp;
        a[((k + 3) - 1) + (k - 1) * lda] = zero;
        Rlartg(a[((k + 1) - 1) + (k - 1) * lda], a[((k + 2) - 1) + (k - 1) * lda], c2, s2, temp);
        a[((k + 1) - 1) + (k - 1) * lda] = temp;
        a[((k + 2) - 1) + (k - 1) * lda] = zero;
        //
        // Apply transformations from the left
        //
        Rrot(istopm - k, &a[((k + 2) - 1) + ((k + 1) - 1) * lda], lda, &a[((k + 3) - 1) + ((k + 1) - 1) * lda], lda, c1, s1);
        Rrot(istopm - k, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda, &a[((k + 2) - 1) + ((k + 1) - 1) * lda], lda, c2, s2);
        //
        Rrot(istopm - k, &b[((k + 2) - 1) + ((k + 1) - 1) * ldb], ldb, &b[((k + 3) - 1) + ((k + 1) - 1) * ldb], ldb, c1, s1);
        Rrot(istopm - k, &b[((k + 1) - 1) + ((k + 1) - 1) * ldb], ldb, &b[((k + 2) - 1) + ((k + 1) - 1) * ldb], ldb, c2, s2);
        if (ilq) {
            Rrot(nq, &q[((k + 2 - qstart + 1) - 1) * ldq], 1, &q[((k + 3 - qstart + 1) - 1) * ldq], 1, c1, s1);
            Rrot(nq, &q[((k + 1 - qstart + 1) - 1) * ldq], 1, &q[((k + 2 - qstart + 1) - 1) * ldq], 1, c2, s2);
        }
        //
    }
    //
    // End of Rlaqz2
    //
}

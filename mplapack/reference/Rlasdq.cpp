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

void Rlasdq(const char *uplo, INTEGER const sqre, INTEGER const n, INTEGER const ncvt, INTEGER const nru, INTEGER const ncc, REAL *d, REAL *e, REAL *vt, INTEGER const ldvt, REAL *u, INTEGER const ldu, REAL *c, INTEGER const ldc, REAL *work, INTEGER &info) {
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    INTEGER iuplo = 0;
    if (Mlsame(uplo, "U")) {
        iuplo = 1;
    }
    if (Mlsame(uplo, "L")) {
        iuplo = 2;
    }
    if (iuplo == 0) {
        info = -1;
    } else if ((sqre < 0) || (sqre > 1)) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ncvt < 0) {
        info = -4;
    } else if (nru < 0) {
        info = -5;
    } else if (ncc < 0) {
        info = -6;
    } else if ((ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < max((INTEGER)1, n))) {
        info = -10;
    } else if (ldu < max((INTEGER)1, nru)) {
        info = -12;
    } else if ((ncc == 0 && ldc < 1) || (ncc > 0 && ldc < max((INTEGER)1, n))) {
        info = -14;
    }
    if (info != 0) {
        Mxerbla("Rlasdq", -info);
        return;
    }
    if (n == 0) {
        return;
    }
    //
    //     ROTATE is true if any singular vectors desired, false otherwise
    //
    bool rotate = (ncvt > 0) || (nru > 0) || (ncc > 0);
    INTEGER np1 = n + 1;
    INTEGER sqre1 = sqre;
    //
    //     If matrix non-square upper bidiagonal, rotate to be lower
    //     bidiagonal.  The rotations are on the right.
    //
    INTEGER i = 0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    REAL r = 0.0;
    const REAL zero = 0.0;
    if ((iuplo == 1) && (sqre1 == 1)) {
        for (i = 1; i <= n - 1; i = i + 1) {
            Rlartg(d[i - 1], e[i - 1], cs, sn, r);
            d[i - 1] = r;
            e[i - 1] = sn * d[(i + 1) - 1];
            d[(i + 1) - 1] = cs * d[(i + 1) - 1];
            if (rotate) {
                work[i - 1] = cs;
                work[(n + i) - 1] = sn;
            }
        }
        Rlartg(d[n - 1], e[n - 1], cs, sn, r);
        d[n - 1] = r;
        e[n - 1] = zero;
        if (rotate) {
            work[n - 1] = cs;
            work[(n + n) - 1] = sn;
        }
        iuplo = 2;
        sqre1 = 0;
        //
        //        Update singular vectors if desired.
        //
        if (ncvt > 0) {
            Rlasr("L", "V", "F", np1, ncvt, &work[1 - 1], &work[np1 - 1], vt, ldvt);
        }
    }
    //
    //     If matrix lower bidiagonal, rotate to be upper bidiagonal
    //     by applying Givens rotations on the left.
    //
    if (iuplo == 2) {
        for (i = 1; i <= n - 1; i = i + 1) {
            Rlartg(d[i - 1], e[i - 1], cs, sn, r);
            d[i - 1] = r;
            e[i - 1] = sn * d[(i + 1) - 1];
            d[(i + 1) - 1] = cs * d[(i + 1) - 1];
            if (rotate) {
                work[i - 1] = cs;
                work[(n + i) - 1] = sn;
            }
        }
        //
        //        If matrix (N+1)-by-N lower bidiagonal, one additional
        //        rotation is needed.
        //
        if (sqre1 == 1) {
            Rlartg(d[n - 1], e[n - 1], cs, sn, r);
            d[n - 1] = r;
            if (rotate) {
                work[n - 1] = cs;
                work[(n + n) - 1] = sn;
            }
        }
        //
        //        Update singular vectors if desired.
        //
        if (nru > 0) {
            if (sqre1 == 0) {
                Rlasr("R", "V", "F", nru, n, &work[1 - 1], &work[np1 - 1], u, ldu);
            } else {
                Rlasr("R", "V", "F", nru, np1, &work[1 - 1], &work[np1 - 1], u, ldu);
            }
        }
        if (ncc > 0) {
            if (sqre1 == 0) {
                Rlasr("L", "V", "F", n, ncc, &work[1 - 1], &work[np1 - 1], c, ldc);
            } else {
                Rlasr("L", "V", "F", np1, ncc, &work[1 - 1], &work[np1 - 1], c, ldc);
            }
        }
    }
    //
    //     Call Rbdsqr to compute the SVD of the reduced real
    //     N-by-N upper bidiagonal matrix.
    //
    Rbdsqr("U", n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info);
    //
    //     Sort the singular values into ascending order (insertion sort on
    //     singular values, but only one transposition per singular vector)
    //
    INTEGER isub = 0;
    REAL smin = 0.0;
    INTEGER j = 0;
    for (i = 1; i <= n; i = i + 1) {
        //
        //        Scan for smallest D(I).
        //
        isub = i;
        smin = d[i - 1];
        for (j = i + 1; j <= n; j = j + 1) {
            if (d[j - 1] < smin) {
                isub = j;
                smin = d[j - 1];
            }
        }
        if (isub != i) {
            //
            //           Swap singular values and vectors.
            //
            d[isub - 1] = d[i - 1];
            d[i - 1] = smin;
            if (ncvt > 0) {
                Rswap(ncvt, &vt[(isub - 1)], ldvt, &vt[(i - 1)], ldvt);
            }
            if (nru > 0) {
                Rswap(nru, &u[(isub - 1) * ldu], 1, &u[(i - 1) * ldu], 1);
            }
            if (ncc > 0) {
                Rswap(ncc, &c[(isub - 1)], ldc, &c[(i - 1)], ldc);
            }
        }
    }
    //
    //     End of Rlasdq
    //
}

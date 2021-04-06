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

void Rlasq1(INTEGER const n, REAL *d, REAL *e, REAL *work, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    REAL sigmn = 0.0;
    REAL sigmx = 0.0;
    if (n < 0) {
        info = -1;
        Mxerbla("Rlasq1", -info);
        return;
    } else if (n == 0) {
        return;
    } else if (n == 1) {
        d[1 - 1] = abs(d[1 - 1]);
        return;
    } else if (n == 2) {
        Rlas2(d[1 - 1], e[1 - 1], d[2 - 1], sigmn, sigmx);
        d[1 - 1] = sigmx;
        d[2 - 1] = sigmn;
        return;
    }
    //
    //     Estimate the largest singular value.
    //
    const REAL zero = 0.0;
    sigmx = zero;
    INTEGER i = 0;
    for (i = 1; i <= n - 1; i = i + 1) {
        d[i - 1] = abs(d[i - 1]);
        sigmx = max(sigmx, abs(e[i - 1]));
    }
    d[n - 1] = abs(d[n - 1]);
    //
    //     Early return if SIGMX is zero (matrix is already diagonal).
    //
    INTEGER iinfo = 0;
    if (sigmx == zero) {
        Rlasrt("D", n, d, iinfo);
        return;
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        sigmx = max(sigmx, d[i - 1]);
    }
    //
    //     Copy D and E into WORK (in the Z format) and scale (squaring the
    //     input data makes scaling by a power of the radix poINTEGERless).
    //
    REAL eps = Rlamch("Precision");
    REAL safmin = Rlamch("Safe minimum");
    REAL scale = sqrt(eps / safmin);
    Rcopy(n, d, 1, &work[1 - 1], 2);
    Rcopy(n - 1, e, 1, &work[2 - 1], 2);
    Rlascl("G", 0, 0, sigmx, scale, 2 * n - 1, 1, work, 2 * n - 1, iinfo);
    //
    //     Compute the q's and e's.
    //
    for (i = 1; i <= 2 * n - 1; i = i + 1) {
        work[i - 1] = pow2(work[i - 1]);
    }
    work[(2 * n) - 1] = zero;
    //
    Rlasq2(n, work, info);
    //
    if (info == 0) {
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = sqrt(work[i - 1]);
        }
        Rlascl("G", 0, 0, scale, sigmx, n, 1, d, n, iinfo);
    } else if (info == 2) {
        //
        //     Maximum number of iterations exceeded.  Move data from WORK
        //     into D and E so the calling subroutine can try to finish
        //
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = sqrt(work[(2 * i - 1) - 1]);
            e[i - 1] = sqrt(work[(2 * i) - 1]);
        }
        Rlascl("G", 0, 0, scale, sigmx, n, 1, d, n, iinfo);
        Rlascl("G", 0, 0, scale, sigmx, n, 1, e, n, iinfo);
    }
    //
    //     End of Rlasq1
    //
}

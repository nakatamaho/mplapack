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

void Rdisna(const char *job, INTEGER const &m, INTEGER const &n, REAL *d, REAL *sep, INTEGER &info) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    bool eigen = Mlsame(job, "E");
    bool left = Mlsame(job, "L");
    bool right = Mlsame(job, "R");
    bool sing = left || right;
    INTEGER k = 0;
    if (eigen) {
        k = m;
    } else if (sing) {
        k = min(m, n);
    }
    bool incr = false;
    bool decr = false;
    INTEGER i = 0;
    const REAL zero = 0.0;
    if (!eigen && !sing) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if (k < 0) {
        info = -3;
    } else {
        incr = true;
        decr = true;
        for (i = 1; i <= k - 1; i = i + 1) {
            if (incr) {
                incr = incr && d[i - 1] <= d[(i + 1) - 1];
            }
            if (decr) {
                decr = decr && d[i - 1] >= d[(i + 1) - 1];
            }
        }
        if (sing && k > 0) {
            if (incr) {
                incr = incr && zero <= d[1 - 1];
            }
            if (decr) {
                decr = decr && d[k - 1] >= zero;
            }
        }
        if (!(incr || decr)) {
            info = -4;
        }
    }
    if (info != 0) {
        Mxerbla("Rdisna", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (k == 0) {
        return;
    }
    //
    //     Compute reciprocal condition numbers
    //
    REAL oldgap = 0.0;
    REAL newgap = 0.0;
    if (k == 1) {
        sep[1 - 1] = dlamch("O");
    } else {
        oldgap = abs(d[2 - 1] - d[1 - 1]);
        sep[1 - 1] = oldgap;
        for (i = 2; i <= k - 1; i = i + 1) {
            newgap = abs(d[(i + 1) - 1] - d[i - 1]);
            sep[i - 1] = min(oldgap, newgap);
            oldgap = newgap;
        }
        sep[k - 1] = oldgap;
    }
    if (sing) {
        if ((left && m > n) || (right && m < n)) {
            if (incr) {
                sep[1 - 1] = min(sep[1 - 1], d[1 - 1]);
            }
            if (decr) {
                sep[k - 1] = min(sep[k - 1], d[k - 1]);
            }
        }
    }
    //
    //     Ensure that reciprocal condition numbers are not less than
    //     threshold, in order to limit the size of the error bound
    //
    REAL eps = dlamch("E");
    REAL safmin = dlamch("S");
    REAL anorm = max(abs(d[1 - 1]), abs(d[k - 1]));
    REAL thresh = 0.0;
    if (anorm == zero) {
        thresh = eps;
    } else {
        thresh = max(eps * anorm, safmin);
    }
    for (i = 1; i <= k; i = i + 1) {
        sep[i - 1] = max(sep[i - 1], thresh);
    }
    //
    //     End of Rdisna
    //
}

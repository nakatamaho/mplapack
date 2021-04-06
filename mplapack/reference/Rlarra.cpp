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

void Rlarra(INTEGER const n, REAL *d, REAL *e, REAL *e2, REAL const spltol, REAL const tnrm, INTEGER &nsplit, INTEGER *isplit, INTEGER &info) {
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
    //
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    //     Compute splitting points
    nsplit = 1;
    const REAL zero = 0.0;
    REAL tmp1 = 0.0;
    INTEGER i = 0;
    REAL eabs = 0.0;
    if (spltol < zero) {
        //        Criterion based on absolute off-diagonal value
        tmp1 = abs(spltol) * tnrm;
        for (i = 1; i <= n - 1; i = i + 1) {
            eabs = abs(e[i - 1]);
            if (eabs <= tmp1) {
                e[i - 1] = zero;
                e2[i - 1] = zero;
                isplit[nsplit - 1] = i;
                nsplit++;
            }
        }
    } else {
        //        Criterion that guarantees relative accuracy
        for (i = 1; i <= n - 1; i = i + 1) {
            eabs = abs(e[i - 1]);
            if (eabs <= spltol * sqrt(abs(d[i - 1])) * sqrt(abs(d[(i + 1) - 1]))) {
                e[i - 1] = zero;
                e2[i - 1] = zero;
                isplit[nsplit - 1] = i;
                nsplit++;
            }
        }
    }
    isplit[nsplit - 1] = n;
    //
    //     End of Rlarra
    //
}

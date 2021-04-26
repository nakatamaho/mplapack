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
#include <mplapack_matgen.h>

COMPLEX
Clatm3(INTEGER const m, INTEGER const n, INTEGER const i, INTEGER const j, INTEGER &isub, INTEGER &jsub, INTEGER const kl, INTEGER const ku, INTEGER const idist, INTEGER *iseed, COMPLEX *d, INTEGER const igrade, COMPLEX *dl, COMPLEX *dr, INTEGER const ipvtng, INTEGER *iwork, REAL const sparse) {
    COMPLEX return_value = (0.0, 0.0);
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //
    //     ..
    //
    //     .. Array Arguments ..
    //
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //
    //     ..
    //
    //     .. Local Scalars ..
    //
    //     ..
    //
    //     .. External Functions ..
    //
    //     ..
    //
    //     .. Intrinsic Functions ..
    //
    //     ..
    //
    //-----------------------------------------------------------------------
    //
    //     .. Executable Statements ..
    //
    //     Check for I and J in range
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (i < 1 || i > m || j < 1 || j > n) {
        isub = i;
        jsub = j;
        return_value = czero;
        return return_value;
    }
    //
    //     Compute subscripts depending on IPVTNG
    //
    if (ipvtng == 0) {
        isub = i;
        jsub = j;
    } else if (ipvtng == 1) {
        isub = iwork[i - 1];
        jsub = j;
    } else if (ipvtng == 2) {
        isub = i;
        jsub = iwork[j - 1];
    } else if (ipvtng == 3) {
        isub = iwork[i - 1];
        jsub = iwork[j - 1];
    }
    //
    //     Check for banding
    //
    if (jsub > isub + ku || jsub < isub - kl) {
        return_value = czero;
        return return_value;
    }
    //
    //     Check for sparsity
    //
    const REAL zero = 0.0;
    if (sparse > zero) {
        if (Rlaran(iseed) < sparse) {
            return_value = czero;
            return return_value;
        }
    }
    //
    //     Compute entry and grade it according to IGRADE
    //
    COMPLEX ctemp = 0.0;
    if (i == j) {
        ctemp = d[i - 1];
    } else {
        ctemp = Clarnd(idist, iseed);
    }
    if (igrade == 1) {
        ctemp = ctemp * dl[i - 1];
    } else if (igrade == 2) {
        ctemp = ctemp * dr[j - 1];
    } else if (igrade == 3) {
        ctemp = ctemp * dl[i - 1] * dr[j - 1];
    } else if (igrade == 4 && i != j) {
        ctemp = ctemp * dl[i - 1] / dl[j - 1];
    } else if (igrade == 5) {
        ctemp = ctemp * dl[i - 1] * conj(dl[j - 1]);
    } else if (igrade == 6) {
        ctemp = ctemp * dl[i - 1] * dl[j - 1];
    }
    return_value = ctemp;
    return return_value;
    //
    //     End of Clatm3
    //
}

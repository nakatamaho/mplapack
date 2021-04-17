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

void Claqgb(INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX *ab, INTEGER const ldab, REAL *r, REAL *c, REAL const rowcnd, REAL const colcnd, REAL const amax, char *equed) {
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (m <= 0 || n <= 0) {
        equed = (char *)"N";
        return;
    }
    //
    //     Initialize LARGE and SMALL.
    //
    REAL small = Rlamch("Safe minimum") / Rlamch("Precision");
    const REAL one = 1.0;
    REAL large = one / small;
    //
    const REAL thresh = 0.1e+0;
    INTEGER j = 0;
    REAL cj = 0.0;
    INTEGER i = 0;
    if (rowcnd >= thresh && amax >= small && amax <= large) {
        //
        //        No row scaling
        //
        if (colcnd >= thresh) {
            //
            //           No column scaling
            //
            equed = (char *)"N";
        } else {
            //
            //           Column scaling
            //
            for (j = 1; j <= n; j = j + 1) {
                cj = c[j - 1];
                for (i = max((INTEGER)1, j - ku); i <= min(m, j + kl); i = i + 1) {
                    ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab] = cj * ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab];
                }
            }
            equed = (char *)"C";
        }
    } else if (colcnd >= thresh) {
        //
        //        Row scaling, no column scaling
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = max((INTEGER)1, j - ku); i <= min(m, j + kl); i = i + 1) {
                ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab] = r[i - 1] * ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab];
            }
        }
        equed = (char *)"R";
    } else {
        //
        //        Row and column scaling
        //
        for (j = 1; j <= n; j = j + 1) {
            cj = c[j - 1];
            for (i = max((INTEGER)1, j - ku); i <= min(m, j + kl); i = i + 1) {
                ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab] = cj * r[i - 1] * ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab];
            }
        }
        equed = (char *)"B";
    }
    //
    //     End of Claqgb
    //
}

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

void Rlaqsp(const char *uplo, INTEGER const n, REAL *ap, REAL *s, REAL const scond, REAL const amax, char *equed) {
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        *equed = 'N';
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
    INTEGER jc = 0;
    INTEGER j = 0;
    REAL cj = 0.0;
    INTEGER i = 0;
    if (scond >= thresh && amax >= small && amax <= large) {
        //
        //        No equilibration
        //
        *equed = 'N';
    } else {
        //
        //        Replace A by diag(S) * A * diag(S).
        //
        if (Mlsame(uplo, "U")) {
            //
            //           Upper triangle of A is stored.
            //
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                cj = s[j - 1];
                for (i = 1; i <= j; i = i + 1) {
                    ap[(jc + i - 1) - 1] = cj * s[i - 1] * ap[(jc + i - 1) - 1];
                }
                jc += j;
            }
        } else {
            //
            //           Lower triangle of A is stored.
            //
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                cj = s[j - 1];
                for (i = j; i <= n; i = i + 1) {
                    ap[(jc + i - j) - 1] = cj * s[i - 1] * ap[(jc + i - j) - 1];
                }
                jc += n - j + 1;
            }
        }
        *equed = 'Y';
    }
    //
    //     End of Rlaqsp
    //
}

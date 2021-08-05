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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Rbdt04(const char *uplo, INTEGER const n, REAL *d, REAL *e, REAL *s, INTEGER const ns, REAL *u, INTEGER const ldu, REAL *vt, INTEGER const ldvt, REAL *work, REAL &resid) {
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    // ======================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible.
    //
    const REAL zero = 0.0;
    resid = zero;
    if (n <= 0 || ns <= 0) {
        return;
    }
    //
    REAL eps = Rlamch("Precision");
    //
    //     Compute S - U' * B * V.
    //
    REAL bnorm = zero;
    //
    INTEGER k = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    if (Mlsame(uplo, "U")) {
        //
        //        B is upper bidiagonal.
        //
        k = 0;
        for (i = 1; i <= ns; i = i + 1) {
            for (j = 1; j <= n - 1; j = j + 1) {
                k++;
                work[k - 1] = d[j - 1] * vt[(i - 1) + (j - 1) * ldvt] + e[j - 1] * vt[(i - 1) + ((j + 1) - 1) * ldvt];
            }
            k++;
            work[k - 1] = d[n - 1] * vt[(i - 1) + (n - 1) * ldvt];
        }
        bnorm = abs(d[1 - 1]);
        for (i = 2; i <= n; i = i + 1) {
            bnorm = max(bnorm, REAL(abs(d[i - 1]) + abs(e[(i - 1) - 1])));
        }
    } else {
        //
        //        B is lower bidiagonal.
        //
        k = 0;
        for (i = 1; i <= ns; i = i + 1) {
            k++;
            work[k - 1] = d[1 - 1] * vt[(i - 1)];
            for (j = 1; j <= n - 1; j = j + 1) {
                k++;
                work[k - 1] = e[j - 1] * vt[(i - 1) + (j - 1) * ldvt] + d[(j + 1) - 1] * vt[(i - 1) + ((j + 1) - 1) * ldvt];
            }
        }
        bnorm = abs(d[n - 1]);
        for (i = 1; i <= n - 1; i = i + 1) {
            bnorm = max(bnorm, REAL(abs(d[i - 1]) + abs(e[i - 1])));
        }
    }
    //
    const REAL one = 1.0;
    Rgemm("T", "N", ns, ns, n, -one, u, ldu, &work[1 - 1], n, zero, &work[(1 + n * ns) - 1], ns);
    //
    //     norm(S - U' * B * V)
    //
    k = n * ns;
    for (i = 1; i <= ns; i = i + 1) {
        work[(k + i) - 1] += s[i - 1];
        resid = max({resid, Rasum(ns, &work[(k + 1) - 1], 1)});
        k += ns;
    }
    //
    if (bnorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        if (bnorm >= resid) {
            resid = (resid / bnorm) / (castREAL(n) * eps);
        } else {
            if (bnorm < one) {
                resid = (min(resid, REAL(castREAL(n) * bnorm)) / bnorm) / (castREAL(n) * eps);
            } else {
                resid = min(REAL(resid / bnorm), castREAL(n)) / (castREAL(n) * eps);
            }
        }
    }
    //
    //     End of Rbdt04
    //
}

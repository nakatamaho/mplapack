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

void Rbdt03(const char *uplo, INTEGER const n, INTEGER const kd, REAL *d, REAL *e, REAL *u, INTEGER const ldu, REAL *s, REAL *vt, INTEGER const ldvt, REAL *work, REAL &resid) {
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
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    resid = zero;
    if (n <= 0) {
        return;
    }
    //
    //     Compute B - U * S * V' one column at a time.
    //
    REAL bnorm = zero;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL one = 1.0;
    if (kd >= 1) {
        //
        //        B is bidiagonal.
        //
        if (Mlsame(uplo, "U")) {
            //
            //           B is upper bidiagonal.
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    work[(n + i) - 1] = s[i - 1] * vt[(i - 1) + (j - 1) * ldvt];
                }
                Rgemv("No transpose", n, n, -one, u, ldu, &work[(n + 1) - 1], 1, zero, work, 1);
                work[j - 1] += d[j - 1];
                if (j > 1) {
                    work[(j - 1) - 1] += e[(j - 1) - 1];
                    bnorm = max(bnorm, abs(d[j - 1]) + abs(e[(j - 1) - 1]));
                } else {
                    bnorm = max(bnorm, abs(d[j - 1]));
                }
                resid = max({resid, Rasum(n, work, 1)});
            }
        } else {
            //
            //           B is lower bidiagonal.
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    work[(n + i) - 1] = s[i - 1] * vt[(i - 1) + (j - 1) * ldvt];
                }
                Rgemv("No transpose", n, n, -one, u, ldu, &work[(n + 1) - 1], 1, zero, work, 1);
                work[j - 1] += d[j - 1];
                if (j < n) {
                    work[(j + 1) - 1] += e[j - 1];
                    bnorm = max(bnorm, abs(d[j - 1]) + abs(e[j - 1]));
                } else {
                    bnorm = max(bnorm, abs(d[j - 1]));
                }
                resid = max({resid, Rasum(n, work, 1)});
            }
        }
    } else {
        //
        //        B is diagonal.
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                work[(n + i) - 1] = s[i - 1] * vt[(i - 1) + (j - 1) * ldvt];
            }
            Rgemv("No transpose", n, n, -one, u, ldu, &work[(n + 1) - 1], 1, zero, work, 1);
            work[j - 1] += d[j - 1];
            resid = max({resid, Rasum(n, work, 1)});
        }
        j = iRamax(n, d, 1);
        bnorm = abs(d[j - 1]);
    }
    //
    //     Compute norm(B - U * S * V') / ( n * norm(B) * EPS )
    //
    REAL eps = Rlamch("Precision");
    //
    if (bnorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        if (bnorm >= resid) {
            resid = (resid / bnorm) / (n.real() * eps);
        } else {
            if (bnorm < one) {
                resid = (min(resid, n.real() * bnorm) / bnorm) / (n.real() * eps);
            } else {
                resid = min(resid / bnorm, n.real()) / (n.real() * eps);
            }
        }
    }
    //
    //     End of Rbdt03
    //
}

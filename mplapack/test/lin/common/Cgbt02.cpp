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

// Derived from LAPACK routine ZGBT02.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Cgbt02(fem::str_cref trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const ldx, COMPLEX *b, INTEGER const ldb, REAL *rwork, REAL &resid) {
    COMPLEX zdum = 0.0;
    //
    // Quick return if N = 0 pr NRHS = 0
    //
    const REAL zero = 0.0;
    if (m <= 0 || n <= 0 || nrhs <= 0) {
        resid = zero;
        return;
    }
    //
    // Exit with RESID = 1/EPS if ANORM = 0.
    //
    REAL eps = Rlamch("Epsilon");
    REAL anorm = zero;
    INTEGER kd = 0;
    INTEGER j = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    REAL temp = 0.0;
    if (Mlsame(trans.elems(), "N")) {
        //
        // Find norm1(A).
        //
        kd = ku + 1;
        for (j = 1; j <= n; j = j + 1) {
            i1 = max(kd + 1 - j, (INTEGER)1);
            i2 = min(kd + m - j, kl + kd);
            if (i2 >= i1) {
                temp = RCasum(i2 - i1 + 1, &a[(i1 - 1) + (j - 1) * lda], 1);
                if (anorm < temp || Risnan(temp)) {
                    anorm = temp;
                }
            }
        }
    } else {
        //
        // Find normI(A).
        //
        for (i1 = 1; i1 <= m; i1 = i1 + 1) {
            rwork[i1 - 1] = zero;
        }
        for (j = 1; j <= n; j = j + 1) {
            kd = ku + 1 - j;
            for (i1 = max((INTEGER)1, j - ku); i1 <= min(m, j + kl); i1 = i1 + 1) {
                rwork[i1 - 1] += cabs1(a[((kd + i1) - 1) + (j - 1) * lda]);
            }
        }
        for (i1 = 1; i1 <= m; i1 = i1 + 1) {
            temp = rwork[i1 - 1];
            if (anorm < temp || Risnan(temp)) {
                anorm = temp;
            }
        }
    }
    const REAL one = 1.0;
    if (anorm <= zero) {
        resid = one / eps;
        return;
    }
    //
    INTEGER n1 = 0;
    if (Mlsame(trans.elems(), "T") || Mlsame(trans.elems(), "C")) {
        n1 = n;
    } else {
        n1 = m;
    }
    //
    // Compute B - op(A)*X
    //
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    for (j = 1; j <= nrhs; j = j + 1) {
        Cgbmv(trans.elems(), m, n, kl, ku, -cone, a, lda, &x[(j - 1) * ldx], 1, cone, &b[(j - 1) * ldb], 1);
    }
    //
    // Compute the maximum over the number of right hand sides of
    // norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).
    //
    resid = zero;
    REAL bnorm = 0.0;
    REAL xnorm = 0.0;
    for (j = 1; j <= nrhs; j = j + 1) {
        bnorm = RCasum(n1, &b[(j - 1) * ldb], 1);
        xnorm = RCasum(n1, &x[(j - 1) * ldx], 1);
        if (xnorm <= zero) {
            resid = one / eps;
        } else {
            resid = max(resid, ((bnorm / anorm) / xnorm) / eps);
        }
    }
    //
    // End of Cgbt02
    //
}

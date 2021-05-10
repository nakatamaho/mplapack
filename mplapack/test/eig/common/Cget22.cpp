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

void Cget22(const char *transa, const char *transe, const char *transw, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *e, INTEGER const lde, COMPLEX *w, COMPLEX *work, REAL *rwork, REAL *result) {
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
    //  =====================================================================
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
    //     Initialize RESULT (in case N=0)
    //
    const REAL zero = 0.0;
    result[1 - 1] = zero;
    result[2 - 1] = zero;
    if (n <= 0) {
        return;
    }
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Precision");
    //
    INTEGER itrnse = 0;
    INTEGER itrnsw = 0;
    char norma;
    char norme;
    //
    if (Mlsame(transa, "T") || Mlsame(transa, "C")) {
        norma = 'I';
    }
    //
    if (Mlsame(transe, "T")) {
        itrnse = 1;
        norme = 'I';
    } else if (Mlsame(transe, "C")) {
        itrnse = 2;
        norme = 'I';
    }
    //
    if (Mlsame(transw, "C")) {
        itrnsw = 1;
    }
    //
    //     Normalization of E:
    //
    const REAL one = 1.0;
    REAL enrmin = one / ulp;
    REAL enrmax = zero;
    INTEGER jvec = 0;
    REAL temp1 = 0.0;
    INTEGER j = 0;
    if (itrnse == 0) {
        for (jvec = 1; jvec <= n; jvec = jvec + 1) {
            temp1 = zero;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = max(temp1, abs(e[(j - 1) + (jvec - 1) * lde].real()) + abs(e[(j - 1) + (jvec - 1) * lde].imag()));
            }
            enrmin = min(enrmin, temp1);
            enrmax = max(enrmax, temp1);
        }
    } else {
        for (jvec = 1; jvec <= n; jvec = jvec + 1) {
            rwork[jvec - 1] = zero;
        }
        //
        for (j = 1; j <= n; j = j + 1) {
            for (jvec = 1; jvec <= n; jvec = jvec + 1) {
                rwork[jvec - 1] = max(rwork[jvec - 1], abs(e[(jvec - 1) + (j - 1) * lde].real()) + abs(e[(jvec - 1) + (j - 1) * lde].imag()));
            }
        }
        //
        for (jvec = 1; jvec <= n; jvec = jvec + 1) {
            enrmin = min(enrmin, rwork[jvec - 1]);
            enrmax = max(enrmax, rwork[jvec - 1]);
        }
    }
    //
    //     Norm of A:
    //
    REAL anorm = max({Clange(&norma, n, n, a, lda, rwork), unfl});
    //
    //     Norm of E:
    //
    REAL enorm = max({Clange(&norme, n, n, e, lde, rwork), ulp});
    //
    //     Norm of error:
    //
    //     Error =  AE - EW
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    Claset("Full", n, n, czero, czero, work, n);
    //
    INTEGER joff = 0;
    INTEGER jcol = 0;
    COMPLEX wtemp = 0.0;
    INTEGER jrow = 0;
    for (jcol = 1; jcol <= n; jcol = jcol + 1) {
        if (itrnsw == 0) {
            wtemp = w[jcol - 1];
        } else {
            wtemp = conj(w[jcol - 1]);
        }
        //
        if (itrnse == 0) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                work[(joff + jrow) - 1] = e[(jrow - 1) + (jcol - 1) * lde] * wtemp;
            }
        } else if (itrnse == 1) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                work[(joff + jrow) - 1] = e[(jcol - 1) + (jrow - 1) * lde] * wtemp;
            }
        } else {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                work[(joff + jrow) - 1] = conj(e[(jcol - 1) + (jrow - 1) * lde]) * wtemp;
            }
        }
        joff += n;
    }
    //
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    Cgemm(transa, transe, n, n, n, cone, a, lda, e, lde, -cone, work, n);
    //
    REAL errnrm = Clange("One", n, n, work, n, rwork) / enorm;
    //
    //     Compute RESULT(1) (avoiding under/overflow)
    //
    if (anorm > errnrm) {
        result[1 - 1] = (errnrm / anorm) / ulp;
    } else {
        if (anorm < one) {
            result[1 - 1] = (min(errnrm, anorm) / anorm) / ulp;
        } else {
            result[1 - 1] = min(errnrm / anorm, one) / ulp;
        }
    }
    //
    //     Compute RESULT(2) : the normalization error in E.
    //
    result[2 - 1] = max(abs(enrmax - one), abs(enrmin - one)) / (castREAL(n) * ulp);
    //
    //     End of Cget22
    //
}

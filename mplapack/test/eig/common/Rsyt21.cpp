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

void Rsyt21(INTEGER const itype, const char *uplo, INTEGER const n, INTEGER const kband, REAL *a, INTEGER const lda, REAL *d, REAL *e, REAL *u, INTEGER const ldu, REAL *v, INTEGER const ldv, REAL *tau, REAL *work, REAL *result) {
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
    const REAL zero = 0.0;
    result[1 - 1] = zero;
    if (itype == 1) {
        result[2 - 1] = zero;
    }
    if (n <= 0) {
        return;
    }
    //
    bool lower = false;
    char cuplo;
    if (Mlsame(uplo, "U")) {
        lower = false;
        cuplo = 'U';
    } else {
        lower = true;
        cuplo = 'L';
    }
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    //
    //     Some Error Checks
    //
    const REAL ten = 10.0;
    if (itype < 1 || itype > 3) {
        result[1 - 1] = ten / ulp;
        return;
    }
    //
    //     Do Test 1
    //
    //     Norm of A:
    //
    const REAL one = 1.0;
    REAL anorm = 0.0;
    if (itype == 3) {
        anorm = one;
    } else {
        anorm = max({Rlansy("1", &cuplo, n, a, lda, work), unfl});
    }
    //
    //     Compute error matrix:
    //
    INTEGER j = 0;
    REAL wnorm = 0.0;
    INTEGER jr = 0;
    REAL vsave = 0.0;
    INTEGER jcol = 0;
    INTEGER jrow = 0;
    INTEGER iinfo = 0;
    if (itype == 1) {
        //
        //        ITYPE=1: error = A - U S U**T
        //
        Rlaset("Full", n, n, zero, zero, work, n);
        Rlacpy(&cuplo, n, n, a, lda, work, n);
        //
        for (j = 1; j <= n; j = j + 1) {
            Rsyr(&cuplo, n, -d[j - 1], &u[(j - 1) * ldu], 1, work, n);
        }
        //
        if (n > 1 && kband == 1) {
            for (j = 1; j <= n - 1; j = j + 1) {
                Rsyr2(&cuplo, n, -e[j - 1], &u[(j - 1) * ldu], 1, &u[((j + 1) - 1) * ldu], 1, work, n);
            }
        }
        wnorm = Rlansy("1", &cuplo, n, work, n, &work[(pow2(n) + 1) - 1]);
        //
    } else if (itype == 2) {
        //
        //        ITYPE=2: error = V S V**T - A
        //
        Rlaset("Full", n, n, zero, zero, work, n);
        //
        if (lower) {
            work[pow2(n) - 1] = d[n - 1];
            for (j = n - 1; j >= 1; j = j - 1) {
                if (kband == 1) {
                    work[((n + 1) * (j - 1) + 2) - 1] = (one - tau[j - 1]) * e[j - 1];
                    for (jr = j + 2; jr <= n; jr = jr + 1) {
                        work[((j - 1) * n + jr) - 1] = -tau[j - 1] * e[j - 1] * v[(jr - 1) + (j - 1) * ldv];
                    }
                }
                //
                vsave = v[((j + 1) - 1) + (j - 1) * ldv];
                v[((j + 1) - 1) + (j - 1) * ldv] = one;
                Rlarfy("L", n - j, &v[((j + 1) - 1) + (j - 1) * ldv], 1, tau[j - 1], &work[((n + 1) * j + 1) - 1], n, &work[(pow2(n) + 1) - 1]);
                v[((j + 1) - 1) + (j - 1) * ldv] = vsave;
                work[((n + 1) * (j - 1) + 1) - 1] = d[j - 1];
            }
        } else {
            work[1 - 1] = d[1 - 1];
            for (j = 1; j <= n - 1; j = j + 1) {
                if (kband == 1) {
                    work[((n + 1) * j) - 1] = (one - tau[j - 1]) * e[j - 1];
                    for (jr = 1; jr <= j - 1; jr = jr + 1) {
                        work[(j * n + jr) - 1] = -tau[j - 1] * e[j - 1] * v[(jr - 1) + ((j + 1) - 1) * ldv];
                    }
                }
                //
                vsave = v[(j - 1) + ((j + 1) - 1) * ldv];
                v[(j - 1) + ((j + 1) - 1) * ldv] = one;
                Rlarfy("U", j, &v[((j + 1) - 1) * ldv], 1, tau[j - 1], work, n, &work[(pow2(n) + 1) - 1]);
                v[(j - 1) + ((j + 1) - 1) * ldv] = vsave;
                work[((n + 1) * j + 1) - 1] = d[(j + 1) - 1];
            }
        }
        //
        for (jcol = 1; jcol <= n; jcol = jcol + 1) {
            if (lower) {
                for (jrow = jcol; jrow <= n; jrow = jrow + 1) {
                    work[(jrow + n * (jcol - 1)) - 1] = work[(jrow + n * (jcol - 1)) - 1] - a[(jrow - 1) + (jcol - 1) * lda];
                }
            } else {
                for (jrow = 1; jrow <= jcol; jrow = jrow + 1) {
                    work[(jrow + n * (jcol - 1)) - 1] = work[(jrow + n * (jcol - 1)) - 1] - a[(jrow - 1) + (jcol - 1) * lda];
                }
            }
        }
        wnorm = Rlansy("1", &cuplo, n, work, n, &work[(pow2(n) + 1) - 1]);
        //
    } else if (itype == 3) {
        //
        //        ITYPE=3: error = U V**T - I
        //
        if (n < 2) {
            return;
        }
        Rlacpy(" ", n, n, u, ldu, work, n);
        if (lower) {
            Rorm2r("R", "T", n, n - 1, n - 1, &v[(2 - 1)], ldv, tau, &work[(n + 1) - 1], n, &work[(pow2(n) + 1) - 1], iinfo);
        } else {
            Rorm2l("R", "T", n, n - 1, n - 1, &v[(2 - 1) * ldv], ldv, tau, work, n, &work[(pow2(n) + 1) - 1], iinfo);
        }
        if (iinfo != 0) {
            result[1 - 1] = ten / ulp;
            return;
        }
        //
        for (j = 1; j <= n; j = j + 1) {
            work[((n + 1) * (j - 1) + 1) - 1] = work[((n + 1) * (j - 1) + 1) - 1] - one;
        }
        //
        wnorm = Rlange("1", n, n, work, n, &work[(pow2(n) + 1) - 1]);
    }
    //
    if (anorm > wnorm) {
        result[1 - 1] = (wnorm / anorm) / (n * ulp);
    } else {
        if (anorm < one) {
            result[1 - 1] = (min(wnorm, n * anorm) / anorm) / (n * ulp);
        } else {
            result[1 - 1] = min(wnorm / anorm, castREAL(n)) / (n * ulp);
        }
    }
    //
    //     Do Test 2
    //
    //     Compute  U U**T - I
    //
    if (itype == 1) {
        Rgemm("N", "C", n, n, n, one, u, ldu, u, ldu, zero, work, n);
        //
        for (j = 1; j <= n; j = j + 1) {
            work[((n + 1) * (j - 1) + 1) - 1] = work[((n + 1) * (j - 1) + 1) - 1] - one;
        }
        //
        result[2 - 1] = min({Rlange("1", n, n, work, n, &work[(pow2(n) + 1) - 1]), castREAL(n)}) / (n * ulp);
    }
    //
    //     End of Rsyt21
    //
}

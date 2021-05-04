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

void Rsbt21(const char *uplo, INTEGER const n, INTEGER const ka, INTEGER const ks, REAL *a, INTEGER const lda, REAL *d, REAL *e, REAL *u, INTEGER const ldu, REAL *work, REAL *result) {
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
    //     Constants
    //
    const REAL zero = 0.0;
    result[1 - 1] = zero;
    result[2 - 1] = zero;
    if (n <= 0) {
        return;
    }
    //
    INTEGER ika = max({(INTEGER)0, min(n - 1, ka)});
    INTEGER lw = (n * (n + 1)) / 2;
    //
    bool lower = false;
    char cuplo[1];
    if (Mlsame(uplo, "U")) {
        lower = false;
        cuplo = "U";
    } else {
        lower = true;
        cuplo = "L";
    }
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    //
    //     Some Error Checks
    //
    //     Do Test 1
    //
    //     Norm of A:
    //
    REAL anorm = max({Rlansb("1", cuplo, n, ika, a, lda, work), unfl});
    //
    //     Compute error matrix:    Error = A - U S U**T
    //
    //     Copy A from SB to SP storage format.
    //
    INTEGER j = 0;
    INTEGER jc = 0;
    INTEGER jr = 0;
    for (jc = 1; jc <= n; jc = jc + 1) {
        if (lower) {
            for (jr = 1; jr <= min(ika + 1, n + 1 - jc); jr = jr + 1) {
                j++;
                work[j - 1] = a[(jr - 1) + (jc - 1) * lda];
            }
            for (jr = ika + 2; jr <= n + 1 - jc; jr = jr + 1) {
                j++;
                work[j - 1] = zero;
            }
        } else {
            for (jr = ika + 2; jr <= jc; jr = jr + 1) {
                j++;
                work[j - 1] = zero;
            }
            for (jr = min(ika, jc - 1); jr >= 0; jr = jr - 1) {
                j++;
                work[j - 1] = a[((ika + 1 - jr) - 1) + (jc - 1) * lda];
            }
        }
    }
    //
    for (j = 1; j <= n; j = j + 1) {
        Rspr(cuplo, n, -d[j - 1], &u[(j - 1) * ldu], 1, work);
    }
    //
    if (n > 1 && ks == 1) {
        for (j = 1; j <= n - 1; j = j + 1) {
            Rspr2(cuplo, n, -e[j - 1], &u[(j - 1) * ldu], 1, &u[((j + 1) - 1) * ldu], 1, work);
        }
    }
    REAL wnorm = Rlansp("1", cuplo, n, work, &work[(lw + 1) - 1]);
    //
    const REAL one = 1.0;
    if (anorm > wnorm) {
        result[1 - 1] = (wnorm / anorm) / (n * ulp);
    } else {
        if (anorm < one) {
            result[1 - 1] = (min(wnorm, n * anorm) / anorm) / (n * ulp);
        } else {
            result[1 - 1] = min(wnorm / anorm, n.real()) / (n * ulp);
        }
    }
    //
    //     Do Test 2
    //
    //     Compute  U U**T - I
    //
    Rgemm("N", "C", n, n, n, one, u, ldu, u, ldu, zero, work, n);
    //
    for (j = 1; j <= n; j = j + 1) {
        work[((n + 1) * (j - 1) + 1) - 1] = work[((n + 1) * (j - 1) + 1) - 1] - one;
    }
    //
    result[2 - 1] = min({Rlange("1", n, n, work, n, &work[(pow2(n) + 1) - 1]), n.real()}) / (n * ulp);
    //
    //     End of Rsbt21
    //
}

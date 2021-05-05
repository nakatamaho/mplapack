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

void Cget54(INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *s, INTEGER const lds, COMPLEX *t, INTEGER const ldt, COMPLEX *u, INTEGER const ldu, COMPLEX *v, INTEGER const ldv, COMPLEX *work, REAL &result) {
    a([lda * star]);
    b([ldb * star]);
    s([lds * star]);
    t([ldt * star]);
    u([ldu * star]);
    v([ldv * star]);
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
    //     .. Local Arrays ..
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
    result = zero;
    if (n <= 0) {
        return;
    }
    //
    //     Constants
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    //
    //     compute the norm of (A,B)
    //
    Clacpy("Full", n, n, a, lda, work, n);
    Clacpy("Full", n, n, b, ldb, &work[(n * n + 1) - 1], n);
    REAL dum[1];
    REAL abnorm = max({Clange("1", n, 2 * n, work, n, dum), unfl});
    //
    //     Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)
    //
    Clacpy(" ", n, n, a, lda, work, n);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    Cgemm("N", "N", n, n, n, cone, u, ldu, s, lds, czero, &work[(n * n + 1) - 1], n);
    //
    Cgemm("N", "C", n, n, n, -cone, &work[(n * n + 1) - 1], n, v, ldv, cone, work, n);
    //
    //     Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)
    //
    Clacpy(" ", n, n, b, ldb, &work[(n * n + 1) - 1], n);
    Cgemm("N", "N", n, n, n, cone, u, ldu, t, ldt, czero, &work[(2 * n * n + 1) - 1], n);
    //
    Cgemm("N", "C", n, n, n, -cone, &work[(2 * n * n + 1) - 1], n, v, ldv, cone, &work[(n * n + 1) - 1], n);
    //
    //     Compute norm(W)/ ( ulp*norm((A,B)) )
    //
    REAL wnorm = Clange("1", n, 2 * n, work, n, dum);
    //
    const REAL one = 1.0;
    if (abnorm > wnorm) {
        result = (wnorm / abnorm) / (2 * n * ulp);
    } else {
        if (abnorm < one) {
            result = (min(wnorm, 2 * n * abnorm) / abnorm) / (2 * n * ulp);
        } else {
            result = min(wnorm / abnorm, (2 * n).real()) / (2 * n * ulp);
        }
    }
    //
    //     End of Cget54
    //
}

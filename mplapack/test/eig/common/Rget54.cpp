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

void Rget54(INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *s, INTEGER const lds, REAL *t, INTEGER const ldt, REAL *u, INTEGER const ldu, REAL *v, INTEGER const ldv, REAL *work, REAL &result) {
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
    Rlacpy("Full", n, n, a, lda, work, n);
    Rlacpy("Full", n, n, b, ldb, &work[(n * n + 1) - 1], n);
    REAL dum[1];
    REAL abnorm = max({Rlange("1", n, 2 * n, work, n, dum), unfl});
    //
    //     Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)
    //
    Rlacpy(" ", n, n, a, lda, work, n);
    const REAL one = 1.0;
    Rgemm("N", "N", n, n, n, one, u, ldu, s, lds, zero, &work[(n * n + 1) - 1], n);
    //
    Rgemm("N", "C", n, n, n, -one, &work[(n * n + 1) - 1], n, v, ldv, one, work, n);
    //
    //     Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)
    //
    Rlacpy(" ", n, n, b, ldb, &work[(n * n + 1) - 1], n);
    Rgemm("N", "N", n, n, n, one, u, ldu, t, ldt, zero, &work[(2 * n * n + 1) - 1], n);
    //
    Rgemm("N", "C", n, n, n, -one, &work[(2 * n * n + 1) - 1], n, v, ldv, one, &work[(n * n + 1) - 1], n);
    //
    //     Compute norm(W)/ ( ulp*norm((A,B)) )
    //
    REAL wnorm = Rlange("1", n, 2 * n, work, n, dum);
    //
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
    //     End of Rget54
    //
}

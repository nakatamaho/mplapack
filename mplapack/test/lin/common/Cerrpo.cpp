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
#include <mplapack_lin.h>

void Cerrpo(const char *path, INTEGER const nunit) {
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER infot;
    INTEGER nout;
    bool ok;
    bool lerr;
    //
    nout = nunit;
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
    INTEGER lda = nmax;
    INTEGER ldaf = nmax;
    COMPLEX b[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX w[2 * nmax];
    COMPLEX x[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
            af[(i - 1) + (j - 1) * ldaf] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
        }
        b[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
    }
    REAL anrm = 1.0;
    ok = true;
    //
    //     Test error exits of the routines that use the Cholesky
    //     decomposition of a Hermitian positive definite matrix.
    //
    INTEGER info = 0;
    REAL r[nmax];
    REAL rcond = 0.0;
    if (Mlsamen(2, c2, "PO")) {
        //
        //        Cpotrf
        //
        infot = 1;
        Cpotrf("/", 0, a, 1, info);
        chkxer("Cpotrf", infot, nout, lerr, ok);
        infot = 2;
        Cpotrf("U", -1, a, 1, info);
        chkxer("Cpotrf", infot, nout, lerr, ok);
        infot = 4;
        Cpotrf("U", 2, a, 1, info);
        chkxer("Cpotrf", infot, nout, lerr, ok);
        //
        //        Cpotf2
        //
        infot = 1;
        Cpotf2("/", 0, a, 1, info);
        chkxer("Cpotf2", infot, nout, lerr, ok);
        infot = 2;
        Cpotf2("U", -1, a, 1, info);
        chkxer("Cpotf2", infot, nout, lerr, ok);
        infot = 4;
        Cpotf2("U", 2, a, 1, info);
        chkxer("Cpotf2", infot, nout, lerr, ok);
        //
        //        Cpotri
        //
        infot = 1;
        Cpotri("/", 0, a, 1, info);
        chkxer("Cpotri", infot, nout, lerr, ok);
        infot = 2;
        Cpotri("U", -1, a, 1, info);
        chkxer("Cpotri", infot, nout, lerr, ok);
        infot = 4;
        Cpotri("U", 2, a, 1, info);
        chkxer("Cpotri", infot, nout, lerr, ok);
        //
        //        Cpotrs
        //
        infot = 1;
        Cpotrs("/", 0, 0, a, 1, b, 1, info);
        chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 2;
        Cpotrs("U", -1, 0, a, 1, b, 1, info);
        chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 3;
        Cpotrs("U", 0, -1, a, 1, b, 1, info);
        chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 5;
        Cpotrs("U", 2, 1, a, 1, b, 2, info);
        chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 7;
        Cpotrs("U", 2, 1, a, 2, b, 1, info);
        chkxer("Cpotrs", infot, nout, lerr, ok);
        //
        //        Cporfs
        //
        infot = 1;
        Cporfs("/", 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 2;
        Cporfs("U", -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 3;
        Cporfs("U", 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 5;
        Cporfs("U", 2, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 7;
        Cporfs("U", 2, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 9;
        Cporfs("U", 2, 1, a, 2, af, 2, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 11;
        Cporfs("U", 2, 1, a, 2, af, 2, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("Cporfs", infot, nout, lerr, ok);
        //
        //        Cpocon
        //
        infot = 1;
        Cpocon("/", 0, a, 1, anrm, rcond, w, r, info);
        chkxer("Cpocon", infot, nout, lerr, ok);
        infot = 2;
        Cpocon("U", -1, a, 1, anrm, rcond, w, r, info);
        chkxer("Cpocon", infot, nout, lerr, ok);
        infot = 4;
        Cpocon("U", 2, a, 1, anrm, rcond, w, r, info);
        chkxer("Cpocon", infot, nout, lerr, ok);
        infot = 5;
        Cpocon("U", 1, a, 1, -anrm, rcond, w, r, info);
        chkxer("Cpocon", infot, nout, lerr, ok);
        //
        //        Cpoequ
        //
        infot = 1;
        Cpoequ(-1, a, 1, r1, rcond, anrm, info);
        chkxer("Cpoequ", infot, nout, lerr, ok);
        infot = 3;
        Cpoequ(2, a, 1, r1, rcond, anrm, info);
        chkxer("Cpoequ", infot, nout, lerr, ok);
        //
        //     Test error exits of the routines that use the Cholesky
        //     decomposition of a Hermitian positive definite packed matrix.
        //
    } else if (Mlsamen(2, c2, "PP")) {
        //
        //        Cpptrf
        //
        infot = 1;
        Cpptrf("/", 0, a, info);
        chkxer("Cpptrf", infot, nout, lerr, ok);
        infot = 2;
        Cpptrf("U", -1, a, info);
        chkxer("Cpptrf", infot, nout, lerr, ok);
        //
        //        Cpptri
        //
        infot = 1;
        Cpptri("/", 0, a, info);
        chkxer("Cpptri", infot, nout, lerr, ok);
        infot = 2;
        Cpptri("U", -1, a, info);
        chkxer("Cpptri", infot, nout, lerr, ok);
        //
        //        Cpptrs
        //
        infot = 1;
        Cpptrs("/", 0, 0, a, b, 1, info);
        chkxer("Cpptrs", infot, nout, lerr, ok);
        infot = 2;
        Cpptrs("U", -1, 0, a, b, 1, info);
        chkxer("Cpptrs", infot, nout, lerr, ok);
        infot = 3;
        Cpptrs("U", 0, -1, a, b, 1, info);
        chkxer("Cpptrs", infot, nout, lerr, ok);
        infot = 6;
        Cpptrs("U", 2, 1, a, b, 1, info);
        chkxer("Cpptrs", infot, nout, lerr, ok);
        //
        //        Cpprfs
        //
        infot = 1;
        Cpprfs("/", 0, 0, a, af, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 2;
        Cpprfs("U", -1, 0, a, af, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 3;
        Cpprfs("U", 0, -1, a, af, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 7;
        Cpprfs("U", 2, 1, a, af, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 9;
        Cpprfs("U", 2, 1, a, af, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("Cpprfs", infot, nout, lerr, ok);
        //
        //        Cppcon
        //
        infot = 1;
        Cppcon("/", 0, a, anrm, rcond, w, r, info);
        chkxer("Cppcon", infot, nout, lerr, ok);
        infot = 2;
        Cppcon("U", -1, a, anrm, rcond, w, r, info);
        chkxer("Cppcon", infot, nout, lerr, ok);
        infot = 4;
        Cppcon("U", 1, a, -anrm, rcond, w, r, info);
        chkxer("Cppcon", infot, nout, lerr, ok);
        //
        //        Cppequ
        //
        infot = 1;
        Cppequ("/", 0, a, r1, rcond, anrm, info);
        chkxer("Cppequ", infot, nout, lerr, ok);
        infot = 2;
        Cppequ("U", -1, a, r1, rcond, anrm, info);
        chkxer("Cppequ", infot, nout, lerr, ok);
        //
        //     Test error exits of the routines that use the Cholesky
        //     decomposition of a Hermitian positive definite band matrix.
        //
    } else if (Mlsamen(2, c2, "PB")) {
        //
        //        Cpbtrf
        //
        infot = 1;
        Cpbtrf("/", 0, 0, a, 1, info);
        chkxer("Cpbtrf", infot, nout, lerr, ok);
        infot = 2;
        Cpbtrf("U", -1, 0, a, 1, info);
        chkxer("Cpbtrf", infot, nout, lerr, ok);
        infot = 3;
        Cpbtrf("U", 1, -1, a, 1, info);
        chkxer("Cpbtrf", infot, nout, lerr, ok);
        infot = 5;
        Cpbtrf("U", 2, 1, a, 1, info);
        chkxer("Cpbtrf", infot, nout, lerr, ok);
        //
        //        Cpbtf2
        //
        infot = 1;
        Cpbtf2("/", 0, 0, a, 1, info);
        chkxer("Cpbtf2", infot, nout, lerr, ok);
        infot = 2;
        Cpbtf2("U", -1, 0, a, 1, info);
        chkxer("Cpbtf2", infot, nout, lerr, ok);
        infot = 3;
        Cpbtf2("U", 1, -1, a, 1, info);
        chkxer("Cpbtf2", infot, nout, lerr, ok);
        infot = 5;
        Cpbtf2("U", 2, 1, a, 1, info);
        chkxer("Cpbtf2", infot, nout, lerr, ok);
        //
        //        Cpbtrs
        //
        infot = 1;
        Cpbtrs("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 2;
        Cpbtrs("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 3;
        Cpbtrs("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 4;
        Cpbtrs("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 6;
        Cpbtrs("U", 2, 1, 1, a, 1, b, 1, info);
        chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 8;
        Cpbtrs("U", 2, 0, 1, a, 1, b, 1, info);
        chkxer("Cpbtrs", infot, nout, lerr, ok);
        //
        //        Cpbrfs
        //
        infot = 1;
        Cpbrfs("/", 0, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 2;
        Cpbrfs("U", -1, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 3;
        Cpbrfs("U", 1, -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 4;
        Cpbrfs("U", 0, 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 6;
        Cpbrfs("U", 2, 1, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 8;
        Cpbrfs("U", 2, 1, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 10;
        Cpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 12;
        Cpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("Cpbrfs", infot, nout, lerr, ok);
        //
        //        Cpbcon
        //
        infot = 1;
        Cpbcon("/", 0, 0, a, 1, anrm, rcond, w, r, info);
        chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 2;
        Cpbcon("U", -1, 0, a, 1, anrm, rcond, w, r, info);
        chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 3;
        Cpbcon("U", 1, -1, a, 1, anrm, rcond, w, r, info);
        chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 5;
        Cpbcon("U", 2, 1, a, 1, anrm, rcond, w, r, info);
        chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 6;
        Cpbcon("U", 1, 0, a, 1, -anrm, rcond, w, r, info);
        chkxer("Cpbcon", infot, nout, lerr, ok);
        //
        //        Cpbequ
        //
        infot = 1;
        Cpbequ("/", 0, 0, a, 1, r1, rcond, anrm, info);
        chkxer("Cpbequ", infot, nout, lerr, ok);
        infot = 2;
        Cpbequ("U", -1, 0, a, 1, r1, rcond, anrm, info);
        chkxer("Cpbequ", infot, nout, lerr, ok);
        infot = 3;
        Cpbequ("U", 1, -1, a, 1, r1, rcond, anrm, info);
        chkxer("Cpbequ", infot, nout, lerr, ok);
        infot = 5;
        Cpbequ("U", 2, 1, a, 1, r1, rcond, anrm, info);
        chkxer("Cpbequ", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrpo
    //
}

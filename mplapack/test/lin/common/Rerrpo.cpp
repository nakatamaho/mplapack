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

void Rerrpo(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
    //
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
    nout = nunit;
    write(nout, star);
    char c2[2] = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    REAL af[nmax * nmax];
    REAL b[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    REAL w[3 * nmax];
    REAL x[nmax];
    INTEGER iw[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / castREAL(i + j);
        }
        b[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        iw[j - 1] = j;
    }
    ok = true;
    //
    INTEGER info = 0;
    REAL anrm = 0.0;
    REAL rcond = 0.0;
    if (Mlsamen(2, c2, "PO")) {
        //
        //        Test error exits of the routines that use the Cholesky
        //        decomposition of a symmetric positive definite matrix.
        //
        //        Rpotrf
        //
        infot = 1;
        Rpotrf("/", 0, a, 1, info);
        chkxer("Rpotrf", infot, nout, lerr, ok);
        infot = 2;
        Rpotrf("U", -1, a, 1, info);
        chkxer("Rpotrf", infot, nout, lerr, ok);
        infot = 4;
        Rpotrf("U", 2, a, 1, info);
        chkxer("Rpotrf", infot, nout, lerr, ok);
        //
        //        Rpotf2
        //
        infot = 1;
        Rpotf2("/", 0, a, 1, info);
        chkxer("Rpotf2", infot, nout, lerr, ok);
        infot = 2;
        Rpotf2("U", -1, a, 1, info);
        chkxer("Rpotf2", infot, nout, lerr, ok);
        infot = 4;
        Rpotf2("U", 2, a, 1, info);
        chkxer("Rpotf2", infot, nout, lerr, ok);
        //
        //        Rpotri
        //
        infot = 1;
        Rpotri("/", 0, a, 1, info);
        chkxer("Rpotri", infot, nout, lerr, ok);
        infot = 2;
        Rpotri("U", -1, a, 1, info);
        chkxer("Rpotri", infot, nout, lerr, ok);
        infot = 4;
        Rpotri("U", 2, a, 1, info);
        chkxer("Rpotri", infot, nout, lerr, ok);
        //
        //        Rpotrs
        //
        infot = 1;
        Rpotrs("/", 0, 0, a, 1, b, 1, info);
        chkxer("Rpotrs", infot, nout, lerr, ok);
        infot = 2;
        Rpotrs("U", -1, 0, a, 1, b, 1, info);
        chkxer("Rpotrs", infot, nout, lerr, ok);
        infot = 3;
        Rpotrs("U", 0, -1, a, 1, b, 1, info);
        chkxer("Rpotrs", infot, nout, lerr, ok);
        infot = 5;
        Rpotrs("U", 2, 1, a, 1, b, 2, info);
        chkxer("Rpotrs", infot, nout, lerr, ok);
        infot = 7;
        Rpotrs("U", 2, 1, a, 2, b, 1, info);
        chkxer("Rpotrs", infot, nout, lerr, ok);
        //
        //        Rporfs
        //
        infot = 1;
        Rporfs("/", 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rporfs", infot, nout, lerr, ok);
        infot = 2;
        Rporfs("U", -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rporfs", infot, nout, lerr, ok);
        infot = 3;
        Rporfs("U", 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rporfs", infot, nout, lerr, ok);
        infot = 5;
        Rporfs("U", 2, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rporfs", infot, nout, lerr, ok);
        infot = 7;
        Rporfs("U", 2, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rporfs", infot, nout, lerr, ok);
        infot = 9;
        Rporfs("U", 2, 1, a, 2, af, 2, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rporfs", infot, nout, lerr, ok);
        infot = 11;
        Rporfs("U", 2, 1, a, 2, af, 2, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rporfs", infot, nout, lerr, ok);
        //
        //        Rpocon
        //
        infot = 1;
        Rpocon("/", 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rpocon", infot, nout, lerr, ok);
        infot = 2;
        Rpocon("U", -1, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rpocon", infot, nout, lerr, ok);
        infot = 4;
        Rpocon("U", 2, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rpocon", infot, nout, lerr, ok);
        //
        //        Rpoequ
        //
        infot = 1;
        Rpoequ(-1, a, 1, r1, rcond, anrm, info);
        chkxer("Rpoequ", infot, nout, lerr, ok);
        infot = 3;
        Rpoequ(2, a, 1, r1, rcond, anrm, info);
        chkxer("Rpoequ", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PP")) {
        //
        //        Test error exits of the routines that use the Cholesky
        //        decomposition of a symmetric positive definite packed matrix.
        //
        //        Rpptrf
        //
        infot = 1;
        Rpptrf("/", 0, a, info);
        chkxer("Rpptrf", infot, nout, lerr, ok);
        infot = 2;
        Rpptrf("U", -1, a, info);
        chkxer("Rpptrf", infot, nout, lerr, ok);
        //
        //        Rpptri
        //
        infot = 1;
        Rpptri("/", 0, a, info);
        chkxer("Rpptri", infot, nout, lerr, ok);
        infot = 2;
        Rpptri("U", -1, a, info);
        chkxer("Rpptri", infot, nout, lerr, ok);
        //
        //        Rpptrs
        //
        infot = 1;
        Rpptrs("/", 0, 0, a, b, 1, info);
        chkxer("Rpptrs", infot, nout, lerr, ok);
        infot = 2;
        Rpptrs("U", -1, 0, a, b, 1, info);
        chkxer("Rpptrs", infot, nout, lerr, ok);
        infot = 3;
        Rpptrs("U", 0, -1, a, b, 1, info);
        chkxer("Rpptrs", infot, nout, lerr, ok);
        infot = 6;
        Rpptrs("U", 2, 1, a, b, 1, info);
        chkxer("Rpptrs", infot, nout, lerr, ok);
        //
        //        Rpprfs
        //
        infot = 1;
        Rpprfs("/", 0, 0, a, af, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rpprfs", infot, nout, lerr, ok);
        infot = 2;
        Rpprfs("U", -1, 0, a, af, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rpprfs", infot, nout, lerr, ok);
        infot = 3;
        Rpprfs("U", 0, -1, a, af, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rpprfs", infot, nout, lerr, ok);
        infot = 7;
        Rpprfs("U", 2, 1, a, af, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rpprfs", infot, nout, lerr, ok);
        infot = 9;
        Rpprfs("U", 2, 1, a, af, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rpprfs", infot, nout, lerr, ok);
        //
        //        Rppcon
        //
        infot = 1;
        Rppcon("/", 0, a, anrm, rcond, w, iw, info);
        chkxer("Rppcon", infot, nout, lerr, ok);
        infot = 2;
        Rppcon("U", -1, a, anrm, rcond, w, iw, info);
        chkxer("Rppcon", infot, nout, lerr, ok);
        //
        //        Rppequ
        //
        infot = 1;
        Rppequ("/", 0, a, r1, rcond, anrm, info);
        chkxer("Rppequ", infot, nout, lerr, ok);
        infot = 2;
        Rppequ("U", -1, a, r1, rcond, anrm, info);
        chkxer("Rppequ", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PB")) {
        //
        //        Test error exits of the routines that use the Cholesky
        //        decomposition of a symmetric positive definite band matrix.
        //
        //        Rpbtrf
        //
        infot = 1;
        Rpbtrf("/", 0, 0, a, 1, info);
        chkxer("Rpbtrf", infot, nout, lerr, ok);
        infot = 2;
        Rpbtrf("U", -1, 0, a, 1, info);
        chkxer("Rpbtrf", infot, nout, lerr, ok);
        infot = 3;
        Rpbtrf("U", 1, -1, a, 1, info);
        chkxer("Rpbtrf", infot, nout, lerr, ok);
        infot = 5;
        Rpbtrf("U", 2, 1, a, 1, info);
        chkxer("Rpbtrf", infot, nout, lerr, ok);
        //
        //        Rpbtf2
        //
        infot = 1;
        Rpbtf2("/", 0, 0, a, 1, info);
        chkxer("Rpbtf2", infot, nout, lerr, ok);
        infot = 2;
        Rpbtf2("U", -1, 0, a, 1, info);
        chkxer("Rpbtf2", infot, nout, lerr, ok);
        infot = 3;
        Rpbtf2("U", 1, -1, a, 1, info);
        chkxer("Rpbtf2", infot, nout, lerr, ok);
        infot = 5;
        Rpbtf2("U", 2, 1, a, 1, info);
        chkxer("Rpbtf2", infot, nout, lerr, ok);
        //
        //        Rpbtrs
        //
        infot = 1;
        Rpbtrs("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("Rpbtrs", infot, nout, lerr, ok);
        infot = 2;
        Rpbtrs("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("Rpbtrs", infot, nout, lerr, ok);
        infot = 3;
        Rpbtrs("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("Rpbtrs", infot, nout, lerr, ok);
        infot = 4;
        Rpbtrs("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("Rpbtrs", infot, nout, lerr, ok);
        infot = 6;
        Rpbtrs("U", 2, 1, 1, a, 1, b, 1, info);
        chkxer("Rpbtrs", infot, nout, lerr, ok);
        infot = 8;
        Rpbtrs("U", 2, 0, 1, a, 1, b, 1, info);
        chkxer("Rpbtrs", infot, nout, lerr, ok);
        //
        //        Rpbrfs
        //
        infot = 1;
        Rpbrfs("/", 0, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        infot = 2;
        Rpbrfs("U", -1, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        infot = 3;
        Rpbrfs("U", 1, -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        infot = 4;
        Rpbrfs("U", 0, 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        infot = 6;
        Rpbrfs("U", 2, 1, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        infot = 8;
        Rpbrfs("U", 2, 1, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        infot = 10;
        Rpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        infot = 12;
        Rpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rpbrfs", infot, nout, lerr, ok);
        //
        //        Rpbcon
        //
        infot = 1;
        Rpbcon("/", 0, 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rpbcon", infot, nout, lerr, ok);
        infot = 2;
        Rpbcon("U", -1, 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rpbcon", infot, nout, lerr, ok);
        infot = 3;
        Rpbcon("U", 1, -1, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rpbcon", infot, nout, lerr, ok);
        infot = 5;
        Rpbcon("U", 2, 1, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rpbcon", infot, nout, lerr, ok);
        //
        //        Rpbequ
        //
        infot = 1;
        Rpbequ("/", 0, 0, a, 1, r1, rcond, anrm, info);
        chkxer("Rpbequ", infot, nout, lerr, ok);
        infot = 2;
        Rpbequ("U", -1, 0, a, 1, r1, rcond, anrm, info);
        chkxer("Rpbequ", infot, nout, lerr, ok);
        infot = 3;
        Rpbequ("U", 1, -1, a, 1, r1, rcond, anrm, info);
        chkxer("Rpbequ", infot, nout, lerr, ok);
        infot = 5;
        Rpbequ("U", 2, 1, a, 1, r1, rcond, anrm, info);
        chkxer("Rpbequ", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrpo
    //
}

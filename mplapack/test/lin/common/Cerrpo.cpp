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

// Derived from LAPACK routine ZERRPO.
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

void Cerrpo(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    write(nout, star);
    fem::str<2> c2 = path(2, 3);
    //
    // Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
    COMPLEX b[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX w[2 * nmax];
    COMPLEX x[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * nmax] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
            af[(i - 1) + (j - 1) * nmax] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
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
    // Test error exits of the routines that use the Cholesky
    // decomposition of a Hermitian positive definite matrix.
    //
    INTEGER info = 0;
    REAL r[nmax];
    REAL rcond = 0.0;
    if (Mlsamen(2, c2.elems, "PO")) {
        //
        // Cpotrf
        //
        srnamt = "Cpotrf";
        infot = 1;
        Cpotrf("/", 0, a, 1, info);
        Chkxer("Cpotrf", infot, nout, lerr, ok);
        infot = 2;
        Cpotrf("U", -1, a, 1, info);
        Chkxer("Cpotrf", infot, nout, lerr, ok);
        infot = 4;
        Cpotrf("U", 2, a, 1, info);
        Chkxer("Cpotrf", infot, nout, lerr, ok);
        //
        // Cpotf2
        //
        srnamt = "Cpotf2";
        infot = 1;
        Cpotf2("/", 0, a, 1, info);
        Chkxer("Cpotf2", infot, nout, lerr, ok);
        infot = 2;
        Cpotf2("U", -1, a, 1, info);
        Chkxer("Cpotf2", infot, nout, lerr, ok);
        infot = 4;
        Cpotf2("U", 2, a, 1, info);
        Chkxer("Cpotf2", infot, nout, lerr, ok);
        //
        // Cpotri
        //
        srnamt = "Cpotri";
        infot = 1;
        Cpotri("/", 0, a, 1, info);
        Chkxer("Cpotri", infot, nout, lerr, ok);
        infot = 2;
        Cpotri("U", -1, a, 1, info);
        Chkxer("Cpotri", infot, nout, lerr, ok);
        infot = 4;
        Cpotri("U", 2, a, 1, info);
        Chkxer("Cpotri", infot, nout, lerr, ok);
        //
        // Cpotrs
        //
        srnamt = "Cpotrs";
        infot = 1;
        Cpotrs("/", 0, 0, a, 1, b, 1, info);
        Chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 2;
        Cpotrs("U", -1, 0, a, 1, b, 1, info);
        Chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 3;
        Cpotrs("U", 0, -1, a, 1, b, 1, info);
        Chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 5;
        Cpotrs("U", 2, 1, a, 1, b, 2, info);
        Chkxer("Cpotrs", infot, nout, lerr, ok);
        infot = 7;
        Cpotrs("U", 2, 1, a, 2, b, 1, info);
        Chkxer("Cpotrs", infot, nout, lerr, ok);
        //
        // Cporfs
        //
        srnamt = "Cporfs";
        infot = 1;
        Cporfs("/", 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 2;
        Cporfs("U", -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 3;
        Cporfs("U", 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 5;
        Cporfs("U", 2, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, r, info);
        Chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 7;
        Cporfs("U", 2, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, r, info);
        Chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 9;
        Cporfs("U", 2, 1, a, 2, af, 2, b, 1, x, 2, r1, r2, w, r, info);
        Chkxer("Cporfs", infot, nout, lerr, ok);
        infot = 11;
        Cporfs("U", 2, 1, a, 2, af, 2, b, 2, x, 1, r1, r2, w, r, info);
        Chkxer("Cporfs", infot, nout, lerr, ok);
        //
        // Cpocon
        //
        srnamt = "Cpocon";
        infot = 1;
        Cpocon("/", 0, a, 1, anrm, rcond, w, r, info);
        Chkxer("Cpocon", infot, nout, lerr, ok);
        infot = 2;
        Cpocon("U", -1, a, 1, anrm, rcond, w, r, info);
        Chkxer("Cpocon", infot, nout, lerr, ok);
        infot = 4;
        Cpocon("U", 2, a, 1, anrm, rcond, w, r, info);
        Chkxer("Cpocon", infot, nout, lerr, ok);
        infot = 5;
        Cpocon("U", 1, a, 1, -anrm, rcond, w, r, info);
        Chkxer("Cpocon", infot, nout, lerr, ok);
        //
        // Cpoequ
        //
        srnamt = "Cpoequ";
        infot = 1;
        Cpoequ(-1, a, 1, r1, rcond, anrm, info);
        Chkxer("Cpoequ", infot, nout, lerr, ok);
        infot = 3;
        Cpoequ(2, a, 1, r1, rcond, anrm, info);
        Chkxer("Cpoequ", infot, nout, lerr, ok);
        //
        // Test error exits of the routines that use the Cholesky
        // decomposition of a Hermitian positive definite packed matrix.
        //
    } else if (Mlsamen(2, c2.elems, "PP")) {
        //
        // Cpptrf
        //
        srnamt = "Cpptrf";
        infot = 1;
        Cpptrf("/", 0, a, info);
        Chkxer("Cpptrf", infot, nout, lerr, ok);
        infot = 2;
        Cpptrf("U", -1, a, info);
        Chkxer("Cpptrf", infot, nout, lerr, ok);
        //
        // Cpptri
        //
        srnamt = "Cpptri";
        infot = 1;
        Cpptri("/", 0, a, info);
        Chkxer("Cpptri", infot, nout, lerr, ok);
        infot = 2;
        Cpptri("U", -1, a, info);
        Chkxer("Cpptri", infot, nout, lerr, ok);
        //
        // Cpptrs
        //
        srnamt = "Cpptrs";
        infot = 1;
        Cpptrs("/", 0, 0, a, b, 1, info);
        Chkxer("Cpptrs", infot, nout, lerr, ok);
        infot = 2;
        Cpptrs("U", -1, 0, a, b, 1, info);
        Chkxer("Cpptrs", infot, nout, lerr, ok);
        infot = 3;
        Cpptrs("U", 0, -1, a, b, 1, info);
        Chkxer("Cpptrs", infot, nout, lerr, ok);
        infot = 6;
        Cpptrs("U", 2, 1, a, b, 1, info);
        Chkxer("Cpptrs", infot, nout, lerr, ok);
        //
        // Cpprfs
        //
        srnamt = "Cpprfs";
        infot = 1;
        Cpprfs("/", 0, 0, a, af, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 2;
        Cpprfs("U", -1, 0, a, af, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 3;
        Cpprfs("U", 0, -1, a, af, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 7;
        Cpprfs("U", 2, 1, a, af, b, 1, x, 2, r1, r2, w, r, info);
        Chkxer("Cpprfs", infot, nout, lerr, ok);
        infot = 9;
        Cpprfs("U", 2, 1, a, af, b, 2, x, 1, r1, r2, w, r, info);
        Chkxer("Cpprfs", infot, nout, lerr, ok);
        //
        // Cppcon
        //
        srnamt = "Cppcon";
        infot = 1;
        Cppcon("/", 0, a, anrm, rcond, w, r, info);
        Chkxer("Cppcon", infot, nout, lerr, ok);
        infot = 2;
        Cppcon("U", -1, a, anrm, rcond, w, r, info);
        Chkxer("Cppcon", infot, nout, lerr, ok);
        infot = 4;
        Cppcon("U", 1, a, -anrm, rcond, w, r, info);
        Chkxer("Cppcon", infot, nout, lerr, ok);
        //
        // Cppequ
        //
        srnamt = "Cppequ";
        infot = 1;
        Cppequ("/", 0, a, r1, rcond, anrm, info);
        Chkxer("Cppequ", infot, nout, lerr, ok);
        infot = 2;
        Cppequ("U", -1, a, r1, rcond, anrm, info);
        Chkxer("Cppequ", infot, nout, lerr, ok);
        //
        // Test error exits of the routines that use the Cholesky
        // decomposition of a Hermitian positive definite band matrix.
        //
    } else if (Mlsamen(2, c2.elems, "PB")) {
        //
        // Cpbtrf
        //
        srnamt = "Cpbtrf";
        infot = 1;
        Cpbtrf("/", 0, 0, a, 1, info);
        Chkxer("Cpbtrf", infot, nout, lerr, ok);
        infot = 2;
        Cpbtrf("U", -1, 0, a, 1, info);
        Chkxer("Cpbtrf", infot, nout, lerr, ok);
        infot = 3;
        Cpbtrf("U", 1, -1, a, 1, info);
        Chkxer("Cpbtrf", infot, nout, lerr, ok);
        infot = 5;
        Cpbtrf("U", 2, 1, a, 1, info);
        Chkxer("Cpbtrf", infot, nout, lerr, ok);
        //
        // Cpbtf2
        //
        srnamt = "Cpbtf2";
        infot = 1;
        Cpbtf2("/", 0, 0, a, 1, info);
        Chkxer("Cpbtf2", infot, nout, lerr, ok);
        infot = 2;
        Cpbtf2("U", -1, 0, a, 1, info);
        Chkxer("Cpbtf2", infot, nout, lerr, ok);
        infot = 3;
        Cpbtf2("U", 1, -1, a, 1, info);
        Chkxer("Cpbtf2", infot, nout, lerr, ok);
        infot = 5;
        Cpbtf2("U", 2, 1, a, 1, info);
        Chkxer("Cpbtf2", infot, nout, lerr, ok);
        //
        // Cpbtrs
        //
        srnamt = "Cpbtrs";
        infot = 1;
        Cpbtrs("/", 0, 0, 0, a, 1, b, 1, info);
        Chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 2;
        Cpbtrs("U", -1, 0, 0, a, 1, b, 1, info);
        Chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 3;
        Cpbtrs("U", 1, -1, 0, a, 1, b, 1, info);
        Chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 4;
        Cpbtrs("U", 0, 0, -1, a, 1, b, 1, info);
        Chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 6;
        Cpbtrs("U", 2, 1, 1, a, 1, b, 1, info);
        Chkxer("Cpbtrs", infot, nout, lerr, ok);
        infot = 8;
        Cpbtrs("U", 2, 0, 1, a, 1, b, 1, info);
        Chkxer("Cpbtrs", infot, nout, lerr, ok);
        //
        // Cpbrfs
        //
        srnamt = "Cpbrfs";
        infot = 1;
        Cpbrfs("/", 0, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 2;
        Cpbrfs("U", -1, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 3;
        Cpbrfs("U", 1, -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 4;
        Cpbrfs("U", 0, 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 6;
        Cpbrfs("U", 2, 1, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 8;
        Cpbrfs("U", 2, 1, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 10;
        Cpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 1, x, 2, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        infot = 12;
        Cpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 2, x, 1, r1, r2, w, r, info);
        Chkxer("Cpbrfs", infot, nout, lerr, ok);
        //
        // Cpbcon
        //
        srnamt = "Cpbcon";
        infot = 1;
        Cpbcon("/", 0, 0, a, 1, anrm, rcond, w, r, info);
        Chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 2;
        Cpbcon("U", -1, 0, a, 1, anrm, rcond, w, r, info);
        Chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 3;
        Cpbcon("U", 1, -1, a, 1, anrm, rcond, w, r, info);
        Chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 5;
        Cpbcon("U", 2, 1, a, 1, anrm, rcond, w, r, info);
        Chkxer("Cpbcon", infot, nout, lerr, ok);
        infot = 6;
        Cpbcon("U", 1, 0, a, 1, -anrm, rcond, w, r, info);
        Chkxer("Cpbcon", infot, nout, lerr, ok);
        //
        // Cpbequ
        //
        srnamt = "Cpbequ";
        infot = 1;
        Cpbequ("/", 0, 0, a, 1, r1, rcond, anrm, info);
        Chkxer("Cpbequ", infot, nout, lerr, ok);
        infot = 2;
        Cpbequ("U", -1, 0, a, 1, r1, rcond, anrm, info);
        Chkxer("Cpbequ", infot, nout, lerr, ok);
        infot = 3;
        Cpbequ("U", 1, -1, a, 1, r1, rcond, anrm, info);
        Chkxer("Cpbequ", infot, nout, lerr, ok);
        infot = 5;
        Cpbequ("U", 2, 1, a, 1, r1, rcond, anrm, info);
        Chkxer("Cpbequ", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Cerrpo
    //
}

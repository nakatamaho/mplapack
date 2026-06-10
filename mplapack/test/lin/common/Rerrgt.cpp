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

// Derived from LAPACK routine DERRGT.
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

void Rerrgt(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    write(nout, star);
    fem::str<2> c2 = path(2, 3);
    const INTEGER nmax = 2;
    REAL d[nmax];
    d[1 - 1] = 1.0;
    d[2 - 1] = 2.0;
    REAL df[nmax];
    df[1 - 1] = 1.0;
    df[2 - 1] = 2.0;
    REAL e[nmax];
    e[1 - 1] = 3.0;
    e[2 - 1] = 4.0;
    REAL ef[nmax];
    ef[1 - 1] = 3.0;
    ef[2 - 1] = 4.0;
    REAL anorm = 1.0;
    ok = true;
    //
    REAL c[nmax];
    REAL f[nmax];
    INTEGER ip[nmax];
    INTEGER info = 0;
    REAL x[nmax];
    REAL cf[nmax];
    REAL b[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    REAL w[nmax];
    INTEGER iw[nmax];
    REAL rcond = 0.0;
    if (Mlsamen(2, c2.elems, "GT")) {
        //
        // Test error exits for the general tridiagonal routines.
        //
        // Rgttrf
        //
        srnamt = "Rgttrf";
        infot = 1;
        Rgttrf(-1, c, d, e, f, ip, info);
        Chkxer("Rgttrf", infot, nout, lerr, ok);
        //
        // Rgttrs
        //
        srnamt = "Rgttrs";
        infot = 1;
        Rgttrs("/", 0, 0, c, d, e, f, ip, x, 1, info);
        Chkxer("Rgttrs", infot, nout, lerr, ok);
        infot = 2;
        Rgttrs("N", -1, 0, c, d, e, f, ip, x, 1, info);
        Chkxer("Rgttrs", infot, nout, lerr, ok);
        infot = 3;
        Rgttrs("N", 0, -1, c, d, e, f, ip, x, 1, info);
        Chkxer("Rgttrs", infot, nout, lerr, ok);
        infot = 10;
        Rgttrs("N", 2, 1, c, d, e, f, ip, x, 1, info);
        Chkxer("Rgttrs", infot, nout, lerr, ok);
        //
        // Rgtrfs
        //
        srnamt = "Rgtrfs";
        infot = 1;
        Rgtrfs("/", 0, 0, c, d, e, cf, df, ef, f, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("Rgtrfs", infot, nout, lerr, ok);
        infot = 2;
        Rgtrfs("N", -1, 0, c, d, e, cf, df, ef, f, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("Rgtrfs", infot, nout, lerr, ok);
        infot = 3;
        Rgtrfs("N", 0, -1, c, d, e, cf, df, ef, f, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("Rgtrfs", infot, nout, lerr, ok);
        infot = 13;
        Rgtrfs("N", 2, 1, c, d, e, cf, df, ef, f, ip, b, 1, x, 2, r1, r2, w, iw, info);
        Chkxer("Rgtrfs", infot, nout, lerr, ok);
        infot = 15;
        Rgtrfs("N", 2, 1, c, d, e, cf, df, ef, f, ip, b, 2, x, 1, r1, r2, w, iw, info);
        Chkxer("Rgtrfs", infot, nout, lerr, ok);
        //
        // Rgtcon
        //
        srnamt = "Rgtcon";
        infot = 1;
        Rgtcon("/", 0, c, d, e, f, ip, anorm, rcond, w, iw, info);
        Chkxer("Rgtcon", infot, nout, lerr, ok);
        infot = 2;
        Rgtcon("I", -1, c, d, e, f, ip, anorm, rcond, w, iw, info);
        Chkxer("Rgtcon", infot, nout, lerr, ok);
        infot = 8;
        Rgtcon("I", 0, c, d, e, f, ip, -anorm, rcond, w, iw, info);
        Chkxer("Rgtcon", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PT")) {
        //
        // Test error exits for the positive definite tridiagonal
        // routines.
        //
        // Rpttrf
        //
        srnamt = "Rpttrf";
        infot = 1;
        Rpttrf(-1, d, e, info);
        Chkxer("Rpttrf", infot, nout, lerr, ok);
        //
        // Rpttrs
        //
        srnamt = "Rpttrs";
        infot = 1;
        Rpttrs(-1, 0, d, e, x, 1, info);
        Chkxer("Rpttrs", infot, nout, lerr, ok);
        infot = 2;
        Rpttrs(0, -1, d, e, x, 1, info);
        Chkxer("Rpttrs", infot, nout, lerr, ok);
        infot = 6;
        Rpttrs(2, 1, d, e, x, 1, info);
        Chkxer("Rpttrs", infot, nout, lerr, ok);
        //
        // Rptrfs
        //
        srnamt = "Rptrfs";
        infot = 1;
        Rptrfs(-1, 0, d, e, df, ef, b, 1, x, 1, r1, r2, w, info);
        Chkxer("Rptrfs", infot, nout, lerr, ok);
        infot = 2;
        Rptrfs(0, -1, d, e, df, ef, b, 1, x, 1, r1, r2, w, info);
        Chkxer("Rptrfs", infot, nout, lerr, ok);
        infot = 8;
        Rptrfs(2, 1, d, e, df, ef, b, 1, x, 2, r1, r2, w, info);
        Chkxer("Rptrfs", infot, nout, lerr, ok);
        infot = 10;
        Rptrfs(2, 1, d, e, df, ef, b, 2, x, 1, r1, r2, w, info);
        Chkxer("Rptrfs", infot, nout, lerr, ok);
        //
        // Rptcon
        //
        srnamt = "Rptcon";
        infot = 1;
        Rptcon(-1, d, e, anorm, rcond, w, info);
        Chkxer("Rptcon", infot, nout, lerr, ok);
        infot = 4;
        Rptcon(0, d, e, -anorm, rcond, w, info);
        Chkxer("Rptcon", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Rerrgt
    //
}

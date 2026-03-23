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

// Derived from LAPACK routine ZERRLS.
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

void Cerrls(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    fem::str<2> c2 = path(2, 3);
    const INTEGER nmax = 2;
    COMPLEX a[nmax * nmax];
    a[0] = COMPLEX(1.0, 0.0);
    a[(2 - 1) * nmax] = COMPLEX(2.0, 0.0);
    a[(2 - 1) + (2 - 1) * nmax] = COMPLEX(3.0, 0.0);
    a[(2 - 1)] = COMPLEX(4.0, 0.0);
    ok = true;
    write(nout, star);
    //
    // Test error exits for the least squares driver routines.
    //
    COMPLEX b[nmax * nmax];
    COMPLEX w[nmax];
    INTEGER info = 0;
    REAL s[nmax];
    REAL rcond = 0.0;
    INTEGER irnk = 0;
    REAL rw[nmax];
    INTEGER ip[nmax];
    if (Mlsamen(2, c2.elems, "LS")) {
        //
        // Cgels
        //
        srnamt = "Cgels";
        infot = 1;
        Cgels("/", 0, 0, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        infot = 2;
        Cgels("N", -1, 0, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        infot = 3;
        Cgels("N", 0, -1, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        infot = 4;
        Cgels("N", 0, 0, -1, a, 1, b, 1, w, 1, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        infot = 6;
        Cgels("N", 2, 0, 0, a, 1, b, 2, w, 2, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        infot = 8;
        Cgels("N", 2, 0, 0, a, 2, b, 1, w, 2, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        infot = 8;
        Cgels("N", 0, 2, 0, a, 1, b, 1, w, 2, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        infot = 10;
        Cgels("N", 1, 1, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgels", infot, nout, lerr, ok);
        //
        // Cgelst
        //
        srnamt = "Cgelst";
        infot = 1;
        Cgelst("/", 0, 0, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        infot = 2;
        Cgelst("N", -1, 0, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        infot = 3;
        Cgelst("N", 0, -1, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        infot = 4;
        Cgelst("N", 0, 0, -1, a, 1, b, 1, w, 1, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        infot = 6;
        Cgelst("N", 2, 0, 0, a, 1, b, 2, w, 2, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        infot = 8;
        Cgelst("N", 2, 0, 0, a, 2, b, 1, w, 2, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        infot = 8;
        Cgelst("N", 0, 2, 0, a, 1, b, 1, w, 2, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        infot = 10;
        Cgelst("N", 1, 1, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgelst", infot, nout, lerr, ok);
        //
        // Cgetsls
        //
        srnamt = "Cgetsls";
        infot = 1;
        Cgetsls("/", 0, 0, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgetsls", infot, nout, lerr, ok);
        infot = 2;
        Cgetsls("N", -1, 0, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgetsls", infot, nout, lerr, ok);
        infot = 3;
        Cgetsls("N", 0, -1, 0, a, 1, b, 1, w, 1, info);
        Chkxer("Cgetsls", infot, nout, lerr, ok);
        infot = 4;
        Cgetsls("N", 0, 0, -1, a, 1, b, 1, w, 1, info);
        Chkxer("Cgetsls", infot, nout, lerr, ok);
        infot = 6;
        Cgetsls("N", 2, 0, 0, a, 1, b, 2, w, 2, info);
        Chkxer("Cgetsls", infot, nout, lerr, ok);
        infot = 8;
        Cgetsls("N", 2, 0, 0, a, 2, b, 1, w, 2, info);
        Chkxer("Cgetsls", infot, nout, lerr, ok);
        infot = 8;
        Cgetsls("N", 0, 2, 0, a, 1, b, 1, w, 2, info);
        Chkxer("Cgetsls", infot, nout, lerr, ok);
        //
        // Cgelss
        //
        srnamt = "Cgelss";
        infot = 1;
        Cgelss(-1, 0, 0, a, 1, b, 1, s, rcond, irnk, w, 1, rw, info);
        Chkxer("Cgelss", infot, nout, lerr, ok);
        infot = 2;
        Cgelss(0, -1, 0, a, 1, b, 1, s, rcond, irnk, w, 1, rw, info);
        Chkxer("Cgelss", infot, nout, lerr, ok);
        infot = 3;
        Cgelss(0, 0, -1, a, 1, b, 1, s, rcond, irnk, w, 1, rw, info);
        Chkxer("Cgelss", infot, nout, lerr, ok);
        infot = 5;
        Cgelss(2, 0, 0, a, 1, b, 2, s, rcond, irnk, w, 2, rw, info);
        Chkxer("Cgelss", infot, nout, lerr, ok);
        infot = 7;
        Cgelss(2, 0, 0, a, 2, b, 1, s, rcond, irnk, w, 2, rw, info);
        Chkxer("Cgelss", infot, nout, lerr, ok);
        //
        // Cgelsy
        //
        srnamt = "Cgelsy";
        infot = 1;
        Cgelsy(-1, 0, 0, a, 1, b, 1, ip, rcond, irnk, w, 10, rw, info);
        Chkxer("Cgelsy", infot, nout, lerr, ok);
        infot = 2;
        Cgelsy(0, -1, 0, a, 1, b, 1, ip, rcond, irnk, w, 10, rw, info);
        Chkxer("Cgelsy", infot, nout, lerr, ok);
        infot = 3;
        Cgelsy(0, 0, -1, a, 1, b, 1, ip, rcond, irnk, w, 10, rw, info);
        Chkxer("Cgelsy", infot, nout, lerr, ok);
        infot = 5;
        Cgelsy(2, 0, 0, a, 1, b, 2, ip, rcond, irnk, w, 10, rw, info);
        Chkxer("Cgelsy", infot, nout, lerr, ok);
        infot = 7;
        Cgelsy(2, 0, 0, a, 2, b, 1, ip, rcond, irnk, w, 10, rw, info);
        Chkxer("Cgelsy", infot, nout, lerr, ok);
        infot = 12;
        Cgelsy(0, 3, 0, a, 1, b, 3, ip, rcond, irnk, w, 1, rw, info);
        Chkxer("Cgelsy", infot, nout, lerr, ok);
        //
        // Cgelsd
        //
        srnamt = "Cgelsd";
        infot = 1;
        Cgelsd(-1, 0, 0, a, 1, b, 1, s, rcond, irnk, w, 10, rw, ip, info);
        Chkxer("Cgelsd", infot, nout, lerr, ok);
        infot = 2;
        Cgelsd(0, -1, 0, a, 1, b, 1, s, rcond, irnk, w, 10, rw, ip, info);
        Chkxer("Cgelsd", infot, nout, lerr, ok);
        infot = 3;
        Cgelsd(0, 0, -1, a, 1, b, 1, s, rcond, irnk, w, 10, rw, ip, info);
        Chkxer("Cgelsd", infot, nout, lerr, ok);
        infot = 5;
        Cgelsd(2, 0, 0, a, 1, b, 2, s, rcond, irnk, w, 10, rw, ip, info);
        Chkxer("Cgelsd", infot, nout, lerr, ok);
        infot = 7;
        Cgelsd(2, 0, 0, a, 2, b, 1, s, rcond, irnk, w, 10, rw, ip, info);
        Chkxer("Cgelsd", infot, nout, lerr, ok);
        infot = 12;
        Cgelsd(2, 2, 1, a, 2, b, 2, s, rcond, irnk, w, 1, rw, ip, info);
        Chkxer("Cgelsd", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Cerrls
    //
}

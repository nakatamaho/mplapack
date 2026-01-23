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

// Derived from LAPACK routine ZERRTZ.
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

void Cerrtz(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    fem::str<2> c2 = path(2, 3);
    const INTEGER nmax = 2;
    COMPLEX a[nmax * nmax];
    a[0] = COMPLEX(1.0, -1.0);
    a[(2 - 1) * nmax] = COMPLEX(2.0, -2.0);
    a[(2 - 1) + (2 - 1) * nmax] = COMPLEX(3.0, -3.0);
    a[(2 - 1)] = COMPLEX(4.0, -4.0);
    COMPLEX w[nmax];
    w[1 - 1] = COMPLEX(0.0, 0.0);
    w[2 - 1] = COMPLEX(0.0, 0.0);
    ok = true;
    //
    // Test error exits for the trapezoidal routines.
    write(nout, star);
    COMPLEX tau[nmax];
    INTEGER info = 0;
    if (Mlsamen(2, c2.elems, "TZ")) {
        //
        // Ctzrzf
        //
        srnamt = "Ctzrzf";
        infot = 1;
        Ctzrzf(-1, 0, a, 1, tau, w, 1, info);
        Chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 2;
        Ctzrzf(1, 0, a, 1, tau, w, 1, info);
        Chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 4;
        Ctzrzf(2, 2, a, 1, tau, w, 1, info);
        Chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 7;
        Ctzrzf(2, 2, a, 2, tau, w, 0, info);
        Chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 7;
        Ctzrzf(2, 3, a, 2, tau, w, 1, info);
        Chkxer("Ctzrzf", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Cerrtz
    //
}

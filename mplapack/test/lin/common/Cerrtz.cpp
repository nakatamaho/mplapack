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
#include <mplapack_debug.h>

void Cerrtz(const char *path, INTEGER const nunit) {
    //
    nout = nunit;
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];

    const INTEGER nmax = 2;
    COMPLEX a[nmax * nmax];
    INTEGER lda = nmax;
    a[(1 - 1) + (1 - 1) * lda] = COMPLEX(1.e+0, -1.e+0);
    a[(1 - 1) + (2 - 1) * lda] = COMPLEX(2.e+0, -2.e+0);
    a[(2 - 1) + (2 - 1) * lda] = COMPLEX(3.e+0, -3.e+0);
    a[(2 - 1)] = COMPLEX(4.e+0, -4.e+0);
    COMPLEX w[nmax];
    w[1 - 1] = COMPLEX(0.e+0, 0.e+0);
    w[2 - 1] = COMPLEX(0.e+0, 0.e+0);
    ok = true;
    //
    //     Test error exits for the trapezoidal routines.
    COMPLEX tau[nmax];
    INTEGER info = 0;
    if (Mlsamen(2, c2, "TZ")) {
        //
        //        Ctzrzf
        //
        strncpy(srnamt, "Ctzrzf", srnamt_len);
        infot = 1;
        Ctzrzf(-1, 0, a, 1, tau, w, 1, info);
        chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 2;
        Ctzrzf(1, 0, a, 1, tau, w, 1, info);
        chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 4;
        Ctzrzf(2, 2, a, 1, tau, w, 1, info);
        chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 7;
        Ctzrzf(2, 2, a, 2, tau, w, 0, info);
        chkxer("Ctzrzf", infot, nout, lerr, ok);
        infot = 7;
        Ctzrzf(2, 3, a, 2, tau, w, 1, info);
        chkxer("Ctzrzf", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrtz
    //
}

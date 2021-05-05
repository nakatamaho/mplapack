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

void Cerrqp(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
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
    char c2[2] = path[(2 - 1) + (3 - 1) * ldpath];
    const INTEGER nmax = 3;
    INTEGER lw = nmax + 1;
    COMPLEX a[nmax * nmax];
    a[(1 - 1)] = COMPLEX(1.0, -1.0);
    a[(2 - 1) * lda] = COMPLEX(2.0e+0, -2.0e+0);
    a[(2 - 1) + (2 - 1) * lda] = COMPLEX(3.0e+0, -3.0e+0);
    a[(2 - 1)] = COMPLEX(4.0e+0, -4.0e+0);
    ok = true;
    write(nout, star);
    //
    //     Test error exits for QR factorization with pivoting
    //
    INTEGER ip[nmax];
    COMPLEX tau[nmax];
    COMPLEX w[2 * nmax + 3 * nmax];
    REAL rw[2 * nmax];
    INTEGER info = 0;
    if (Mlsamen(2, c2, "QP")) {
        //
        //        Cgeqp3
        //
        infot = 1;
        Cgeqp3(-1, 0, a, 1, ip, tau, w, lw, rw, info);
        chkxer("Cgeqp3", infot, nout, lerr, ok);
        infot = 2;
        Cgeqp3(1, -1, a, 1, ip, tau, w, lw, rw, info);
        chkxer("Cgeqp3", infot, nout, lerr, ok);
        infot = 4;
        Cgeqp3(2, 3, a, 1, ip, tau, w, lw, rw, info);
        chkxer("Cgeqp3", infot, nout, lerr, ok);
        infot = 8;
        Cgeqp3(2, 2, a, 2, ip, tau, w, lw - 10, rw, info);
        chkxer("Cgeqp3", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrqp
    //
}

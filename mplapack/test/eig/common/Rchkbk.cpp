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

void Rchkbk(INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    char buf[1024];
    INTEGER lmax[2];
    INTEGER ninfo = 0;
    INTEGER knt = 0;
    const REAL zero = 0.0;
    REAL rmax = 0.0;
    REAL eps = 0.0;
    REAL safmin = 0.0;
    INTEGER n = 0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    const INTEGER lde = 20;
    REAL scale[lde];
    INTEGER i = 0;
    REAL e[lde * lde];
    INTEGER j = 0;
    REAL ein[lde * lde];
    INTEGER ldein = lde;
    INTEGER info = 0;
    REAL vmax = 0.0;
    REAL x = 0.0;
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    // ======================================================================
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
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    ninfo = 0;
    knt = 0;
    rmax = zero;
    eps = Rlamch("E");
    safmin = Rlamch("S");
    double dtmp;
    //
statement_10:
    //
    read(nin, star), n, ilo, ihi;
    if (n == 0) {
        goto statement_60;
    }
    //
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= n; i = i + 1) {
            rloop, dtmp;
            scale[i - 1] = dtmp;
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, dtmp;
                e[(i - 1) + (j - 1) * lde] = dtmp;
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, dtmp;
                ein[(i - 1) + (j - 1) * ldein] = dtmp;
            }
        }
    }
    //
    knt++;
    Rgebak("B", "R", n, ilo, ihi, scale, n, e, lde, info);
    //
    if (info != 0) {
        ninfo++;
        lmax[1 - 1] = knt;
    }
    //
    vmax = zero;
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            x = abs(e[(i - 1) + (j - 1) * lde] - ein[(i - 1) + (j - 1) * ldein]) / eps;
            if (abs(e[(i - 1) + (j - 1) * lde]) > safmin) {
                x = x / abs(e[(i - 1) + (j - 1) * lde]);
            }
            vmax = max(vmax, x);
        }
    }
    //
    if (vmax > rmax) {
        lmax[2 - 1] = knt;
        rmax = vmax;
    }
    //
    goto statement_10;
//
statement_60:
    //
    write(nout, "(1x,'.. test output of Rgebak .. ')");
    //
    sprintnum_short(buf, rmax);
    write(nout, "(1x,'value of largest test error             = ',a)"), buf;
    write(nout, "(1x,'example number where info is not zero   = ',i4)"), lmax[1 - 1];
    write(nout, "(1x,'example number having largest error     = ',i4)"), lmax[2 - 1];
    write(nout, "(1x,'number of examples where info is not 0  = ',i4)"), ninfo;
    write(nout, "(1x,'total number of examples tested         = ',i4)"), knt;
    //
    //     End of Rchkbk
    //
}

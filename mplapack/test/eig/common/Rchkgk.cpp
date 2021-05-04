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

void Rchkgk(INTEGER const nin, INTEGER const nout) {
    common_read read(cmn);
    common_write write(cmn);
    arr_1d<4, int> lmax;
    INTEGER ninfo = 0;
    INTEGER knt = 0;
    const REAL zero = 0.0;
    REAL rmax = 0.0;
    REAL eps = 0.0;
    INTEGER n = 0;
    INTEGER m = 0;
    INTEGER i = 0;
    const INTEGER lda = 50;
    INTEGER j = 0;
    const INTEGER ldb = 50;
    const INTEGER ldvl = 50;
    const INTEGER ldvr = 50;
    const INTEGER ldwork = 50;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    arr_1d<lda, REAL> lscale;
    arr_1d<lda, REAL> rscale;
    INTEGER info = 0;
    const REAL one = 1.0;
    const INTEGER lde = 50;
    const INTEGER ldf = 50;
    REAL vmax = 0.0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialization
    //
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    lmax[3 - 1] = 0;
    lmax[4 - 1] = 0;
    ninfo = 0;
    knt = 0;
    rmax = zero;
    //
    eps = Rlamch("Precision");
//
statement_10:
    read(nin, star), n, m;
    if (n == 0) {
        goto statement_100;
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, a(i, j);
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, b(i, j);
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= m; j = j + 1) {
                rloop, vl(i, j);
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= m; j = j + 1) {
                rloop, vr(i, j);
            }
        }
    }
    //
    knt++;
    //
    anorm = Rlange("M", n, n, a, lda, work);
    bnorm = Rlange("M", n, n, b, ldb, work);
    //
    Rlacpy("FULL", n, n, a, lda, af, lda);
    Rlacpy("FULL", n, n, b, ldb, bf, ldb);
    //
    Rggbal("B", n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
    if (info != 0) {
        ninfo++;
        lmax[1 - 1] = knt;
    }
    //
    Rlacpy("FULL", n, m, vl, ldvl, vlf, ldvl);
    Rlacpy("FULL", n, m, vr, ldvr, vrf, ldvr);
    //
    Rggbak("B", "L", n, ilo, ihi, lscale, rscale, m, vl, ldvl, info);
    if (info != 0) {
        ninfo++;
        lmax[2 - 1] = knt;
    }
    //
    Rggbak("B", "R", n, ilo, ihi, lscale, rscale, m, vr, ldvr, info);
    if (info != 0) {
        ninfo++;
        lmax[3 - 1] = knt;
    }
    //
    //     Test of Rggbak
    //
    //     Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
    //     where tilde(A) denotes the transformed matrix.
    //
    Rgemm("N", "N", n, m, n, one, af, lda, vr, ldvr, zero, work, ldwork);
    Rgemm("T", "N", m, m, n, one, vl, ldvl, work, ldwork, zero, e, lde);
    //
    Rgemm("N", "N", n, m, n, one, a, lda, vrf, ldvr, zero, work, ldwork);
    Rgemm("T", "N", m, m, n, one, vlf, ldvl, work, ldwork, zero, f, ldf);
    //
    vmax = zero;
    for (j = 1; j <= m; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            vmax = max(vmax, abs(e[(i - 1) + (j - 1) * lde] - f[(i - 1) + (j - 1) * ldf]));
        }
    }
    vmax = vmax / (eps * max(anorm, bnorm));
    if (vmax > rmax) {
        lmax[4 - 1] = knt;
        rmax = vmax;
    }
    //
    //     Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR
    //
    Rgemm("N", "N", n, m, n, one, bf, ldb, vr, ldvr, zero, work, ldwork);
    Rgemm("T", "N", m, m, n, one, vl, ldvl, work, ldwork, zero, e, lde);
    //
    Rgemm("N", "N", n, m, n, one, b, ldb, vrf, ldvr, zero, work, ldwork);
    Rgemm("T", "N", m, m, n, one, vlf, ldvl, work, ldwork, zero, f, ldf);
    //
    vmax = zero;
    for (j = 1; j <= m; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            vmax = max(vmax, abs(e[(i - 1) + (j - 1) * lde] - f[(i - 1) + (j - 1) * ldf]));
        }
    }
    vmax = vmax / (eps * max(anorm, bnorm));
    if (vmax > rmax) {
        lmax[4 - 1] = knt;
        rmax = vmax;
    }
    //
    goto statement_10;
//
statement_100:
    //
    write(nout, "(1x,'.. test output of Rggbak .. ')");
    //
    write(nout, "(' value of largest test error                  =',d12.3)"), rmax;
    write(nout, "(' example number where Rggbal info is not 0    =',i4)"), lmax(1);
    write(nout, "(' example number where Rggbak(L) info is not 0 =',i4)"), lmax(2);
    write(nout, "(' example number where Rggbak(R) info is not 0 =',i4)"), lmax(3);
    write(nout, "(' example number having largest error          =',i4)"), lmax(4);
    write(nout, "(' number of examples where info is not 0       =',i4)"), ninfo;
    write(nout, "(' total number of examples tested              =',i4)"), knt;
    //
    //     End of Rchkgk
    //
}

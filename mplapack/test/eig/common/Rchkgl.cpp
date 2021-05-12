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

void Rchkgl(INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    char buf[1024];
    INTEGER lmax[5];
    INTEGER ninfo = 0;
    INTEGER knt = 0;
    const REAL zero = 0.0;
    REAL rmax = 0.0;
    REAL eps = 0.0;
    INTEGER n = 0;
    INTEGER i = 0;
    const INTEGER lda = 20;
    REAL a[lda * lda];
    INTEGER j = 0;
    const INTEGER ldb = 20;
    REAL b[ldb * ldb];
    INTEGER iloin = 0;
    INTEGER ihiin = 0;
    REAL ain[lda * lda];
    REAL bin[ldb * ldb];
    INTEGER ldain = lda;
    INTEGER ldbin = ldb;
    REAL lsclin[lda];
    REAL rsclin[lda];
    const INTEGER lwork = 6 * lda;
    REAL work[lwork];
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL lscale[lda];
    REAL rscale[lda];
    INTEGER info = 0;
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
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    lmax[3 - 1] = 0;
    ninfo = 0;
    knt = 0;
    rmax = zero;
    double dtmp;
    //
    eps = Rlamch("Precision");
//
statement_10:
    //
    read(nin, star), n;
    if (n == 0) {
        goto statement_90;
    }
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, dtmp;
                a[(i - 1) + (j - 1) * lda] = dtmp;
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, dtmp;
                b[(i - 1) + (j - 1) * ldb] = dtmp;
            }
        }
    }
    //
    read(nin, star), iloin, ihiin;
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, dtmp;
                ain[(i - 1) + (j - 1) * ldain] = dtmp;
            }
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, dtmp;
                bin[(i - 1) + (j - 1) * ldbin] = dtmp;
            }
        }
    }
    //
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= n; i = i + 1) {
            rloop, dtmp;
            lsclin[i - 1] = dtmp;
        }
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= n; i = i + 1) {
            rloop, dtmp;
            rsclin[i - 1] = dtmp;
        }
    }
    //
    anorm = Rlange("M", n, n, a, lda, work);
    bnorm = Rlange("M", n, n, b, ldb, work);
    //
    knt++;
    //
    Rggbal("B", n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
    //
    if (info != 0) {
        ninfo++;
        lmax[1 - 1] = knt;
    }
    //
    if (ilo != iloin || ihi != ihiin) {
        ninfo++;
        lmax[2 - 1] = knt;
    }
    //
    vmax = zero;
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            vmax = max(vmax, abs(a[(i - 1) + (j - 1) * lda] - ain[(i - 1) + (j - 1) * ldain]));
            vmax = max(vmax, abs(b[(i - 1) + (j - 1) * ldb] - bin[(i - 1) + (j - 1) * ldbin]));
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        vmax = max(vmax, abs(lscale[i - 1] - lsclin[i - 1]));
        vmax = max(vmax, abs(rscale[i - 1] - rsclin[i - 1]));
    }
    //
    vmax = vmax / (eps * max(anorm, bnorm));
    //
    if (vmax > rmax) {
        lmax[3 - 1] = knt;
        rmax = vmax;
    }
    //
    goto statement_10;
//
statement_90:
    //
    write(nout, "(1x,'.. test output of Rggbal .. ')");
    //
    sprintnum_short(buf, rmax);
    write(nout, "(1x,'value of largest test error            = ',a)"), buf;
    write(nout, "(1x,'example number where info is not zero  = ',i4)"), lmax[1 - 1];
    write(nout, "(1x,'example number where ILO or IHI wrong  = ',i4)"), lmax[2 - 1];
    write(nout, "(1x,'example number having largest error    = ',i4)"), lmax[3 - 1];
    write(nout, "(1x,'number of examples where info is not 0 = ',i4)"), ninfo;
    write(nout, "(1x,'total number of examples tested        = ',i4)"), knt;
    //
    //     End of Rchkgl
    //
}

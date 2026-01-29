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

// Derived from LAPACK routine ZCHKGK.
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
#include <mplapack_eig.h>
#include <memory>

void Cchkgk(INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    COMPLEX cdum = 0.0;
    INTEGER lmax[4];
    INTEGER ninfo = 0;
    INTEGER knt = 0;
    const REAL zero = 0.0;
    REAL rmax = 0.0;
    REAL eps = 0.0;
    INTEGER n = 0;
    INTEGER m = 0;
    INTEGER i = 0;
    const INTEGER lda = 50;
    auto a_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, lda * lda));
    COMPLEX *a = a_storage.get();
    INTEGER j = 0;
    const INTEGER ldb = 50;
    auto b_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldb * ldb));
    COMPLEX *b = b_storage.get();
    const INTEGER ldvl = 50;
    auto vl_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldvl * ldvl));
    COMPLEX *vl = vl_storage.get();
    const INTEGER ldvr = 50;
    auto vr_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldvr * ldvr));
    COMPLEX *vr = vr_storage.get();
    const INTEGER lrwork = 6 * 50;
    REAL rwork[lrwork];
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    auto af_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, lda * lda));
    COMPLEX *af = af_storage.get();
    auto bf_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldb * ldb));
    COMPLEX *bf = bf_storage.get();
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL lscale[lda];
    REAL rscale[lda];
    INTEGER info = 0;
    auto vlf_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldvl * ldvl));
    COMPLEX *vlf = vlf_storage.get();
    auto vrf_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldvr * ldvr));
    COMPLEX *vrf = vrf_storage.get();
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const INTEGER ldwork = 50;
    auto work_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldwork * ldwork));
    COMPLEX *work = work_storage.get();
    const INTEGER lde = 50;
    auto e_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, lde * lde));
    COMPLEX *e = e_storage.get();
    const INTEGER ldf = 50;
    auto f_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, ldf * ldf));
    COMPLEX *f = f_storage.get();
    REAL vmax = 0.0;
    static const char *format_9999 = "(1x,'.. test output of Cggbak .. ')";
    static const char *format_9998 = "(' value of largest test error                  =',d12.3)";
    static const char *format_9997 = "(' example number where Cggbal info is not 0    =',i4)";
    static const char *format_9996 = "(' example number where Cggbak(L) info is not 0 =',i4)";
    static const char *format_9995 = "(' example number where Cggbak(R) info is not 0 =',i4)";
    static const char *format_9994 = "(' example number having largest error          =',i4)";
    static const char *format_9992 = "(' number of examples where info is not 0       =',i4)";
    static const char *format_9991 = "(' total number of examples tested              =',i4)";
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
                rloop, a[(i - 1) + (j - 1) * lda];
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, b[(i - 1) + (j - 1) * ldb];
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= m; j = j + 1) {
                rloop, vl[(i - 1) + (j - 1) * ldvl];
            }
        }
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= m; j = j + 1) {
                rloop, vr[(i - 1) + (j - 1) * ldvr];
            }
        }
    }
    //
    knt++;
    //
    anorm = Clange("M", n, n, a, lda, rwork);
    bnorm = Clange("M", n, n, b, ldb, rwork);
    //
    Clacpy("FULL", n, n, a, lda, af, lda);
    Clacpy("FULL", n, n, b, ldb, bf, ldb);
    //
    Cggbal("B", n, a, lda, b, ldb, ilo, ihi, lscale, rscale, rwork, info);
    if (info != 0) {
        ninfo++;
        lmax[1 - 1] = knt;
    }
    //
    Clacpy("FULL", n, m, vl, ldvl, vlf, ldvl);
    Clacpy("FULL", n, m, vr, ldvr, vrf, ldvr);
    //
    Cggbak("B", "L", n, ilo, ihi, lscale, rscale, m, vl, ldvl, info);
    if (info != 0) {
        ninfo++;
        lmax[2 - 1] = knt;
    }
    //
    Cggbak("B", "R", n, ilo, ihi, lscale, rscale, m, vr, ldvr, info);
    if (info != 0) {
        ninfo++;
        lmax[3 - 1] = knt;
    }
    //
    // Test of Cggbak
    //
    // Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
    // where tilde(A) denotes the transformed matrix.
    //
    Cgemm("N", "N", n, m, n, cone, af, lda, vr, ldvr, czero, work, ldwork);
    Cgemm("C", "N", m, m, n, cone, vl, ldvl, work, ldwork, czero, e, lde);
    //
    Cgemm("N", "N", n, m, n, cone, a, lda, vrf, ldvr, czero, work, ldwork);
    Cgemm("C", "N", m, m, n, cone, vlf, ldvl, work, ldwork, czero, f, ldf);
    //
    vmax = zero;
    for (j = 1; j <= m; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            vmax = max(vmax, cabs1(e[(i - 1) + (j - 1) * lde] - f[(i - 1) + (j - 1) * ldf]));
        }
    }
    vmax = vmax / (eps * max(anorm, bnorm));
    if (vmax > rmax) {
        lmax[4 - 1] = knt;
        rmax = vmax;
    }
    //
    // Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR
    //
    Cgemm("N", "N", n, m, n, cone, bf, ldb, vr, ldvr, czero, work, ldwork);
    Cgemm("C", "N", m, m, n, cone, vl, ldvl, work, ldwork, czero, e, lde);
    //
    Cgemm("n", "n", n, m, n, cone, b, ldb, vrf, ldvr, czero, work, ldwork);
    Cgemm("C", "N", m, m, n, cone, vlf, ldvl, work, ldwork, czero, f, ldf);
    //
    vmax = zero;
    for (j = 1; j <= m; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            vmax = max(vmax, cabs1(e[(i - 1) + (j - 1) * lde] - f[(i - 1) + (j - 1) * ldf]));
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
    write(nout, format_9999);
    //
    write(nout, format_9998), rmax;
    write(nout, format_9997), lmax[1 - 1];
    write(nout, format_9996), lmax[2 - 1];
    write(nout, format_9995), lmax[3 - 1];
    write(nout, format_9994), lmax[4 - 1];
    write(nout, format_9992), ninfo;
    write(nout, format_9991), knt;
    //
    // End of Cchkgk
    //
}

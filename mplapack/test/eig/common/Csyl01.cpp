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

// Derived from LAPACK routine ZSYL01.
// Original LAPACK authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>
#include <memory>

void Csyl01(REAL const thresh, INTEGER *nfail, REAL *rmax, INTEGER *ninfo, INTEGER &knt) {
    INTEGER allocatestatus = 0;
    const INTEGER maxm = 185;
    std::unique_ptr<COMPLEX[]> a_storage;
    COMPLEX *a = nullptr;
    a_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, maxm * maxm));
    a = a_storage.get();
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    allocatestatus = 0;
    const INTEGER maxn = 192;
    std::unique_ptr<COMPLEX[]> b_storage;
    COMPLEX *b = nullptr;
    b_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, maxn * maxn));
    b = b_storage.get();
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    allocatestatus = 0;
    std::unique_ptr<COMPLEX[]> c_storage;
    COMPLEX *c = nullptr;
    c_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, maxm * maxn));
    c = c_storage.get();
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    allocatestatus = 0;
    std::unique_ptr<COMPLEX[]> cc_storage;
    COMPLEX *cc = nullptr;
    cc_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, maxm * maxn));
    cc = cc_storage.get();
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    allocatestatus = 0;
    std::unique_ptr<COMPLEX[]> x_storage;
    COMPLEX *x = nullptr;
    x_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, maxm * maxn));
    x = x_storage.get();
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    allocatestatus = 0;
    const INTEGER ldswork = 36;
    std::unique_ptr<REAL[]> swork_storage;
    REAL *swork = nullptr;
    swork_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldswork * 103));
    swork = swork_storage.get();
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    //
    // Get machine parameters
    //
    REAL eps = Rlamch("P");
    REAL smlnum = Rlamch("S") / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    // Expect INFO = 0
    REAL vm[2];
    vm[1 - 1] = one;
    // Expect INFO = 1
    vm[2 - 1] = 0.05;
    //
    // Begin test loop
    //
    ninfo[1 - 1] = 0;
    ninfo[2 - 1] = 0;
    nfail[1 - 1] = 0;
    nfail[2 - 1] = 0;
    nfail[3 - 1] = 0;
    const REAL zero = 0.0;
    rmax[1 - 1] = zero;
    rmax[2 - 1] = zero;
    knt = 0;
    INTEGER iseed[4];
    iseed[1 - 1] = 1;
    iseed[2 - 1] = 1;
    iseed[3 - 1] = 1;
    iseed[4 - 1] = 1;
    REAL scale = one;
    REAL scale3 = one;
    INTEGER j = 0;
    INTEGER isgn = 0;
    INTEGER m = 0;
    INTEGER kla = 0;
    INTEGER kua = 0;
    auto d_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, max(maxm, maxn)));
    COMPLEX *d = d_storage.get();
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    COMPLEX duml[maxm];
    COMPLEX dumr[maxn];
    INTEGER iwork[maxm + maxn + 2];
    INTEGER iinfo = 0;
    INTEGER i = 0;
    REAL dum[maxn];
    REAL anrm = 0.0;
    INTEGER n = 0;
    INTEGER klb = 0;
    INTEGER kub = 0;
    REAL bnrm = 0.0;
    REAL tnrm = 0.0;
    INTEGER itrana = 0;
    fem::str<1> trana;
    INTEGER itranb = 0;
    fem::str<1> tranb;
    REAL xnrm = 0.0;
    COMPLEX rmul = 0.0;
    REAL res1 = 0.0;
    REAL res = 0.0;
    INTEGER info = 0;
    for (j = 1; j <= 2; j = j + 1) {
        for (isgn = -1; isgn <= 1; isgn = isgn + 2) {
            // Reset seed (overwritten by LATMR)
            iseed[1 - 1] = 1;
            iseed[2 - 1] = 1;
            iseed[3 - 1] = 1;
            iseed[4 - 1] = 1;
            for (m = 32; m <= maxm; m = m + 51) {
                kla = 0;
                kua = m - 1;
                Clatmr(m, m, "S", iseed, "N", d, 6, one, cone, "T", "N", duml, 1, one, dumr, 1, one, "N", iwork, kla, kua, zero, one, "NO", a, maxm, iwork, iinfo);
                for (i = 1; i <= m; i = i + 1) {
                    a[(i - 1) + (i - 1) * maxm] = a[(i - 1) + (i - 1) * maxm] * vm[j - 1];
                }
                anrm = Clange("M", m, m, a, maxm, dum);
                for (n = 51; n <= maxn; n = n + 47) {
                    klb = 0;
                    kub = n - 1;
                    Clatmr(n, n, "S", iseed, "N", d, 6, one, cone, "T", "N", duml, 1, one, dumr, 1, one, "N", iwork, klb, kub, zero, one, "NO", b, maxn, iwork, iinfo);
                    for (i = 1; i <= n; i = i + 1) {
                        b[(i - 1) + (i - 1) * maxn] = b[(i - 1) + (i - 1) * maxn] * vm[j - 1];
                    }
                    bnrm = Clange("M", n, n, b, maxn, dum);
                    tnrm = max(anrm, bnrm);
                    Clatmr(m, n, "S", iseed, "N", d, 6, one, cone, "T", "N", duml, 1, one, dumr, 1, one, "N", iwork, m, n, zero, one, "NO", c, maxm, iwork, iinfo);
                    for (itrana = 1; itrana <= 2; itrana = itrana + 1) {
                        if (itrana == 1) {
                            trana = "N";
                        }
                        if (itrana == 2) {
                            trana = "C";
                        }
                        for (itranb = 1; itranb <= 2; itranb = itranb + 1) {
                            if (itranb == 1) {
                                tranb = "N";
                            }
                            if (itranb == 2) {
                                tranb = "C";
                            }
                            knt++;
                            //
                            Clacpy("All", m, n, c, maxm, x, maxm);
                            Clacpy("All", m, n, c, maxm, cc, maxm);
                            Ctrsyl(trana.elems, tranb.elems, isgn, m, n, a, maxm, b, maxn, x, maxm, scale, iinfo);
                            if (iinfo != 0) {
                                ninfo[1 - 1]++;
                            }
                            xnrm = Clange("M", m, n, x, maxm, dum);
                            rmul = cone;
                            if (xnrm > one && tnrm > one) {
                                if (xnrm > bignum / tnrm) {
                                    rmul = cone / max(xnrm, tnrm);
                                }
                            }
                            Cgemm(trana.elems, "N", m, n, m, rmul, a, maxm, x, maxm, -scale * rmul, cc, maxm);
                            Cgemm("N", tranb.elems, m, n, n, castREAL(isgn) * rmul, x, maxm, b, maxn, cone, cc, maxm);
                            res1 = Clange("M", m, n, cc, maxm, dum);
                            res = res1 / max(smlnum, smlnum * xnrm, ((abs(rmul) * tnrm) * eps) * xnrm);
                            if (res > thresh) {
                                nfail[1 - 1]++;
                            }
                            if (res > rmax[1 - 1]) {
                                rmax[1 - 1] = res;
                            }
                            //
                            Clacpy("All", m, n, c, maxm, x, maxm);
                            Clacpy("All", m, n, c, maxm, cc, maxm);
                            Ctrsyl3(trana.elems, tranb.elems, isgn, m, n, a, maxm, b, maxn, x, maxm, scale3, swork, ldswork, info);
                            if (info != 0) {
                                ninfo[2 - 1]++;
                            }
                            xnrm = Clange("M", m, n, x, maxm, dum);
                            rmul = cone;
                            if (xnrm > one && tnrm > one) {
                                if (xnrm > bignum / tnrm) {
                                    rmul = cone / max(xnrm, tnrm);
                                }
                            }
                            Cgemm(trana.elems, "N", m, n, m, rmul, a, maxm, x, maxm, -scale3 * rmul, cc, maxm);
                            Cgemm("N", tranb.elems, m, n, n, castREAL(isgn) * rmul, x, maxm, b, maxn, cone, cc, maxm);
                            res1 = Clange("M", m, n, cc, maxm, dum);
                            res = res1 / max(smlnum, smlnum * xnrm, ((abs(rmul) * tnrm) * eps) * xnrm);
                            // Verify that TRSYL3 only flushes if TRSYL flushes (but
                            // there may be cases where TRSYL3 avoid flushing).
                            if (scale3 == zero && scale > zero || iinfo != info) {
                                nfail[3 - 1]++;
                            }
                            if (res > thresh || Risnan(res)) {
                                nfail[2 - 1]++;
                            }
                            if (res > rmax[2 - 1]) {
                                rmax[2 - 1] = res;
                            }
                        }
                    }
                }
            }
        }
    }
    //
    // End of Csyl01
    //
}

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

void Rget35(REAL &rmax, INTEGER &lmax, INTEGER &ninfo, INTEGER &knt) {
    INTEGER *idim(sve.idim, [8]);
    arr_ref<int, 3> ival(sve.ival, [6 * 6 * 8]);
    //
    if (is_called_first_time) {
        {
            static const INTEGER values[] = {1, 2, 3, 4, 3, 3, 6, 4};
            data_of_type<int>(FEM_VALUES_AND_SIZE), idim;
        }
        {
            data_values data;
            data.values, 1, 35 * datum(0), 1, 2, 4 * datum(0), -2, 0, 28 * datum(0);
            data.values, 1, 5 * datum(0), 5, 1, 2, 3 * datum(0), -8, -2;
            data.values, 1, 21 * datum(0), 3, 4, 4 * datum(0), -5, 3, 4 * datum(0);
            data.values, 1, 2, 1, 4, 2 * datum(0), -3, -9, -1;
            data.values, 1, 14 * datum(0), 1, 5 * datum(0), 2, 3, 4 * datum(0), 5;
            data.values, 6, 7, 21 * datum(0), 1, 5 * datum(0), 1, 3, -4;
            data.values, 3 * datum(0), 2, 5, 2, 21 * datum(0), 1, 2, 4 * datum(0);
            data.values, -2, 0, 4 * datum(0), 5, 6, 3, 4, 2 * datum(0);
            data.values, -1, -9, -5, 2, 2 * datum(0), 4 * datum(8), 5, 6;
            data.values, 4 * datum(9), -7, 5, 1, 5 * datum(0), 1, 5, 2;
            data.values, 3 * datum(0), 2, -21, 5, 3 * datum(0), 1, 2, 3;
            data.values, 4, 14 * datum(0);
            data, ival;
        }
    }
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
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Get machine parameters
    //
    REAL eps = Rlamch("P");
    const REAL four = 4.0;
    REAL smlnum = Rlamch("S") * four / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Set up test case parameters
    //
    REAL vm1[3];
    vm1[1 - 1] = sqrt(smlnum);
    vm1[2 - 1] = one;
    vm1[3 - 1] = sqrt(bignum);
    REAL vm2[3];
    vm2[1 - 1] = one;
    const REAL two = 2.0;
    vm2[2 - 1] = one + two * eps;
    vm2[3 - 1] = two;
    //
    knt = 0;
    ninfo = 0;
    lmax = 0;
    const REAL zero = 0.0;
    rmax = zero;
    //
    //     Begin test loop
    //
    INTEGER itrana = 0;
    INTEGER itranb = 0;
    INTEGER isgn = 0;
    INTEGER ima = 0;
    INTEGER imlda1 = 0;
    INTEGER imlda2 = 0;
    INTEGER imloff = 0;
    INTEGER imb = 0;
    INTEGER imldb1 = 0;
    char trana;
    char tranb;
    INTEGER m = 0;
    INTEGER n = 0;
    REAL tnrm = 0.0;
    INTEGER i = 0;
    INTEGER j = 0;
    REAL a[6 * 6];
    REAL b[6 * 6];
    REAL cnrm = 0.0;
    REAL c[6 * 6];
    REAL cc[6 * 6];
    REAL scale = 0.0;
    INTEGER info = 0;
    REAL dum[1];
    REAL xnrm = 0.0;
    REAL rmul = 0.0;
    REAL res1 = 0.0;
    REAL res = 0.0;
    for (itrana = 1; itrana <= 2; itrana = itrana + 1) {
        for (itranb = 1; itranb <= 2; itranb = itranb + 1) {
            for (isgn = -1; isgn <= 1; isgn = isgn + 2) {
                for (ima = 1; ima <= 8; ima = ima + 1) {
                    for (imlda1 = 1; imlda1 <= 3; imlda1 = imlda1 + 1) {
                        for (imlda2 = 1; imlda2 <= 3; imlda2 = imlda2 + 1) {
                            for (imloff = 1; imloff <= 2; imloff = imloff + 1) {
                                for (imb = 1; imb <= 8; imb = imb + 1) {
                                    for (imldb1 = 1; imldb1 <= 3; imldb1 = imldb1 + 1) {
                                        if (itrana == 1) {
                                            trana = "N";
                                        }
                                        if (itrana == 2) {
                                            trana = "T";
                                        }
                                        if (itranb == 1) {
                                            tranb = "N";
                                        }
                                        if (itranb == 2) {
                                            tranb = "T";
                                        }
                                        m = idim[ima - 1];
                                        n = idim[imb - 1];
                                        tnrm = zero;
                                        for (i = 1; i <= m; i = i + 1) {
                                            for (j = 1; j <= m; j = j + 1) {
                                                a[(i - 1) + (j - 1) * lda] = ivali, j, ima;
                                                if (abs(i - j) <= 1) {
                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * vm1[imlda1 - 1];
                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * vm2[imlda2 - 1];
                                                } else {
                                                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * vm1[imloff - 1];
                                                }
                                                tnrm = max(tnrm, abs(a[(i - 1) + (j - 1) * lda]));
                                            }
                                        }
                                        for (i = 1; i <= n; i = i + 1) {
                                            for (j = 1; j <= n; j = j + 1) {
                                                b[(i - 1) + (j - 1) * ldb] = ivali, j, imb;
                                                if (abs(i - j) <= 1) {
                                                    b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] * vm1[imldb1 - 1];
                                                } else {
                                                    b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] * vm1[imloff - 1];
                                                }
                                                tnrm = max(tnrm, abs(b[(i - 1) + (j - 1) * ldb]));
                                            }
                                        }
                                        cnrm = zero;
                                        for (i = 1; i <= m; i = i + 1) {
                                            for (j = 1; j <= n; j = j + 1) {
                                                c[(i - 1) + (j - 1) * ldc] = sin((i * j).real());
                                                cnrm = max(cnrm, &c[(i - 1) + (j - 1) * ldc]);
                                                cc[(i - 1) + (j - 1) * ldcc] = c[(i - 1) + (j - 1) * ldc];
                                            }
                                        }
                                        knt++;
                                        Rtrsyl(trana, tranb, isgn, m, n, a, 6, b, 6, c, 6, scale, info);
                                        if (info != 0) {
                                            ninfo++;
                                        }
                                        xnrm = Rlange("M", m, n, c, 6, dum);
                                        rmul = one;
                                        if (xnrm > one && tnrm > one) {
                                            if (xnrm > bignum / tnrm) {
                                                rmul = one / max(xnrm, tnrm);
                                            }
                                        }
                                        Rgemm(trana, "N", m, n, m, rmul, a, 6, c, 6, -scale * rmul, cc, 6);
                                        Rgemm("N", tranb, m, n, n, isgn.real() * rmul, c, 6, b, 6, one, cc, 6);
                                        res1 = Rlange("M", m, n, cc, 6, dum);
                                        res = res1 / max({smlnum, smlnum * xnrm, ((rmul * tnrm) * eps) * xnrm});
                                        if (res > rmax) {
                                            lmax = knt;
                                            rmax = res;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //
    //     End of Rget35
    //
}

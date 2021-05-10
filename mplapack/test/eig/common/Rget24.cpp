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

INTEGER seldim, selopt;
bool selval[20];
REAL selwi[20], selwr[20];

bool _Rslect(REAL const zr, REAL const zi) {
    bool return_value = false;
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
    //     .. Arrays in Common ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Parameters ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    const REAL zero = 0.0;
    REAL rmin = 0.0;
    INTEGER i = 0;
    REAL x = 0.0;
    if (selopt == 0) {
        return_value = (zr < zero);
    } else {
        rmin = Rlapy2(zr - selwr[1 - 1], zi - selwi[1 - 1]);
        return_value = selval[1 - 1];
        for (i = 2; i <= seldim; i = i + 1) {
            x = Rlapy2(zr - selwr[i - 1], zi - selwi[i - 1]);
            if (x <= rmin) {
                rmin = x;
                return_value = selval[i - 1];
            }
        }
    }
    return return_value;
    //
    //     End of Rslect
    //
}

void Rget24(bool const comp, INTEGER const jtype, REAL const thresh, INTEGER *iseed, INTEGER const nounit, INTEGER const n, REAL *a, INTEGER const lda, REAL *h, REAL *ht, REAL *wr, REAL *wi, REAL *wrt, REAL *wit, REAL *wrtmp, REAL *witmp, REAL *vs, INTEGER const ldvs, REAL *vs1, REAL const rcdein, REAL const rcdvin, INTEGER const nslct, INTEGER *islct, REAL *result, REAL *work, INTEGER const lwork, INTEGER *iwork, bool *bwork, INTEGER &info) {

    INTEGER ldh = lda;
    INTEGER ldht = lda;
    INTEGER ldvs1 = ldvs;
    common cmn;
    common_write write(cmn);
    //
    const REAL zero = 0.0;
    INTEGER i = 0;
    const REAL one = 1.0;
    REAL smlnum = 0.0;
    REAL ulp = 0.0;
    REAL ulpinv = 0.0;
    INTEGER liwork = 0;
    INTEGER isort = 0;
    char sort;
    INTEGER rsub = 0;
    INTEGER sdim = 0;
    REAL rconde = 0.0;
    REAL rcondv = 0.0;
    INTEGER iinfo = 0;
    INTEGER j = 0;
    REAL anorm = 0.0;
    REAL wnorm = 0.0;
    REAL tmp = 0.0;
    INTEGER knteig = 0;
    INTEGER sdim1 = 0;
    REAL rcnde1 = 0.0;
    REAL rcndv1 = 0.0;
    const REAL epsin = 5.9605e-8;
    REAL eps = 0.0;
    INTEGER ipnt[20];
    INTEGER kmin = 0;
    REAL vrmin = 0.0;
    REAL vimin = 0.0;
    INTEGER itmp = 0;
    REAL v = 0.0;
    REAL tol = 0.0;
    REAL tolin = 0.0;
    static const char *format_9998 = "(' Rget24: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',3(i5,','),i5,')')";
    static const char *format_9999 = "(' Rget24: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
                                     "', INPUT EXAMPLE NUMBER = ',i4)";
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
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
    //     .. Arrays in Common ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Check for errors
    //
    info = 0;
    if (thresh < zero) {
        info = -3;
    } else if (nounit <= 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (lda < 1 || lda < n) {
        info = -8;
    } else if (ldvs < 1 || ldvs < n) {
        info = -18;
    } else if (lwork < 3 * n) {
        info = -26;
    }
    //
    if (info != 0) {
        Mxerbla("Rget24", -info);
        return;
    }
    //
    //     Quick return if nothing to do
    //
    for (i = 1; i <= 17; i = i + 1) {
        result[i - 1] = -one;
    }
    //
    if (n == 0) {
        return;
    }
    //
    //     Important constants
    //
    smlnum = Rlamch("Safe minimum");
    ulp = Rlamch("Precision");
    ulpinv = one / ulp;
    //
    //     Perform tests (1)-(13)
    //
    selopt = 0;
    liwork = n * n;
    for (isort = 0; isort <= 1; isort = isort + 1) {
        if (isort == 0) {
            sort = 'N';
            rsub = 0;
        } else {
            sort = 'S';
            rsub = 6;
        }
        //
        //        Compute Schur form and Schur vectors, and test them
        //
        Rlacpy("F", n, n, a, lda, h, lda);
        Rgeesx("V", &sort, _Rslect, "N", n, h, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[(1 + rsub) - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx1", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx1", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            return;
        }
        if (isort == 0) {
            Rcopy(n, wr, 1, wrtmp, 1);
            Rcopy(n, wi, 1, witmp, 1);
        }
        //
        //        Do Test (1) or Test (7)
        //
        result[(1 + rsub) - 1] = zero;
        for (j = 1; j <= n - 2; j = j + 1) {
            for (i = j + 2; i <= n; i = i + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != zero) {
                    result[(1 + rsub) - 1] = ulpinv;
                }
            }
        }
        for (i = 1; i <= n - 2; i = i + 1) {
            if (h[((i + 1) - 1) + (i - 1) * ldh] != zero && h[((i + 2) - 1) + ((i + 1) - 1) * ldh] != zero) {
                result[(1 + rsub) - 1] = ulpinv;
            }
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            if (h[((i + 1) - 1) + (i - 1) * ldh] != zero) {
                if (h[(i - 1) + (i - 1) * ldh] != h[((i + 1) - 1) + ((i + 1) - 1) * ldh] || h[(i - 1) + ((i + 1) - 1) * ldh] == zero || sign(one, h[((i + 1) - 1) + (i - 1) * ldh]) == sign(one, h[(i - 1) + ((i + 1) - 1) * ldh])) {
                    result[(1 + rsub) - 1] = ulpinv;
                }
            }
        }
        //
        //        Test (2) or (8): Compute norm(A - Q*H*Q') / (norm(A) * N * ULP)
        //
        //        Copy A to VS1, used as workspace
        //
        Rlacpy(" ", n, n, a, lda, vs1, ldvs);
        //
        //        Compute Q*H and store in HT.
        //
        Rgemm("No transpose", "No transpose", n, n, n, one, vs, ldvs, h, lda, zero, ht, lda);
        //
        //        Compute A - Q*H*Q'
        //
        Rgemm("No transpose", "Transpose", n, n, n, -one, ht, lda, vs, ldvs, one, vs1, ldvs);
        //
        anorm = max({Rlange("1", n, n, a, lda, work), smlnum});
        wnorm = Rlange("1", n, n, vs1, ldvs, work);
        //
        if (anorm > wnorm) {
            result[(2 + rsub) - 1] = (wnorm / anorm) / (n * ulp);
        } else {
            if (anorm < one) {
                result[(2 + rsub) - 1] = (min(wnorm, n * anorm) / anorm) / (n * ulp);
            } else {
                result[(2 + rsub) - 1] = min(wnorm / anorm, castREAL(n)) / (n * ulp);
            }
        }
        //
        //        Test (3) or (9):  Compute norm( I - Q'*Q ) / ( N * ULP )
        //
        Rort01("Columns", n, n, vs, ldvs, work, lwork, result[(3 + rsub) - 1]);
        //
        //        Do Test (4) or Test (10)
        //
        result[(4 + rsub) - 1] = zero;
        for (i = 1; i <= n; i = i + 1) {
            if (h[(i - 1) + (i - 1) * ldh] != wr[i - 1]) {
                result[(4 + rsub) - 1] = ulpinv;
            }
        }
        if (n > 1) {
            if (h[(2 - 1)] == zero && wi[1 - 1] != zero) {
                result[(4 + rsub) - 1] = ulpinv;
            }
            if (h[(n - 1) + ((n - 1) - 1) * ldh] == zero && wi[n - 1] != zero) {
                result[(4 + rsub) - 1] = ulpinv;
            }
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            if (h[((i + 1) - 1) + (i - 1) * ldh] != zero) {
                tmp = sqrt(abs(h[((i + 1) - 1) + (i - 1) * ldh])) * sqrt(abs(h[(i - 1) + ((i + 1) - 1) * ldh]));
                result[(4 + rsub) - 1] = max({result[(4 + rsub) - 1], abs(wi[i - 1] - tmp) / max(ulp * tmp, smlnum)});
                result[(4 + rsub) - 1] = max({result[(4 + rsub) - 1], abs(wi[(i + 1) - 1] + tmp) / max(ulp * tmp, smlnum)});
            } else if (i > 1) {
                if (h[((i + 1) - 1) + (i - 1) * ldh] == zero && h[(i - 1) + ((i - 1) - 1) * ldh] == zero && wi[i - 1] != zero) {
                    result[(4 + rsub) - 1] = ulpinv;
                }
            }
        }
        //
        //        Do Test (5) or Test (11)
        //
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("N", &sort, _Rslect, "N", n, ht, lda, sdim, wrt, wit, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[(5 + rsub) - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx2", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx2", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            goto statement_250;
        }
        //
        result[(5 + rsub) - 1] = zero;
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != ht[(i - 1) + (j - 1) * ldht]) {
                    result[(5 + rsub) - 1] = ulpinv;
                }
            }
        }
        //
        //        Do Test (6) or Test (12)
        //
        result[(6 + rsub) - 1] = zero;
        for (i = 1; i <= n; i = i + 1) {
            if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                result[(6 + rsub) - 1] = ulpinv;
            }
        }
        //
        //        Do Test (13)
        //
        if (isort == 1) {
            result[13 - 1] = zero;
            knteig = 0;
            for (i = 1; i <= n; i = i + 1) {
                if (_Rslect(wr[i - 1], wi[i - 1]) || _Rslect(wr[i - 1], -wi[i - 1])) {
                    knteig++;
                }
                if (i < n) {
                    if ((_Rslect(wr[(i + 1) - 1], wi[(i + 1) - 1]) || _Rslect(wr[(i + 1) - 1], -wi[(i + 1) - 1])) && (!(_Rslect(wr[i - 1], wi[i - 1]) || _Rslect(wr[i - 1], -wi[i - 1]))) && iinfo != n + 2) {
                        result[13 - 1] = ulpinv;
                    }
                }
            }
            if (sdim != knteig) {
                result[13 - 1] = ulpinv;
            }
        }
        //
    }
    //
    //     If there is enough workspace, perform tests (14) and (15)
    //     as well as (10) through (13)
    //
    if (lwork >= n + (n * n) / 2) {
        //
        //        Compute both RCONDE and RCONDV with VS
        //
        sort = 'S';
        result[14 - 1] = zero;
        result[15 - 1] = zero;
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("V", &sort, _Rslect, "B", n, ht, lda, sdim1, wrt, wit, vs1, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[14 - 1] = ulpinv;
            result[15 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx3", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx3", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            goto statement_250;
        }
        //
        //        Perform tests (10), (11), (12), and (13)
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                result[10 - 1] = ulpinv;
            }
            for (j = 1; j <= n; j = j + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != ht[(i - 1) + (j - 1) * ldht]) {
                    result[11 - 1] = ulpinv;
                }
                if (vs[(i - 1) + (j - 1) * ldvs] != vs1[(i - 1) + (j - 1) * ldvs1]) {
                    result[12 - 1] = ulpinv;
                }
            }
        }
        if (sdim != sdim1) {
            result[13 - 1] = ulpinv;
        }
        //
        //        Compute both RCONDE and RCONDV without VS, and compare
        //
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("N", &sort, _Rslect, "B", n, ht, lda, sdim1, wrt, wit, vs1, ldvs, rcnde1, rcndv1, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[14 - 1] = ulpinv;
            result[15 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx4", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx4", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            goto statement_250;
        }
        //
        //        Perform tests (14) and (15)
        //
        if (rcnde1 != rconde) {
            result[14 - 1] = ulpinv;
        }
        if (rcndv1 != rcondv) {
            result[15 - 1] = ulpinv;
        }
        //
        //        Perform tests (10), (11), (12), and (13)
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                result[10 - 1] = ulpinv;
            }
            for (j = 1; j <= n; j = j + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != ht[(i - 1) + (j - 1) * ldht]) {
                    result[11 - 1] = ulpinv;
                }
                if (vs[(i - 1) + (j - 1) * ldvs] != vs1[(i - 1) + (j - 1) * ldvs1]) {
                    result[12 - 1] = ulpinv;
                }
            }
        }
        if (sdim != sdim1) {
            result[13 - 1] = ulpinv;
        }
        //
        //        Compute RCONDE with VS, and compare
        //
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("V", &sort, _Rslect, "E", n, ht, lda, sdim1, wrt, wit, vs1, ldvs, rcnde1, rcndv1, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[14 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx5", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx5", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            goto statement_250;
        }
        //
        //        Perform test (14)
        //
        if (rcnde1 != rconde) {
            result[14 - 1] = ulpinv;
        }
        //
        //        Perform tests (10), (11), (12), and (13)
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                result[10 - 1] = ulpinv;
            }
            for (j = 1; j <= n; j = j + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != ht[(i - 1) + (j - 1) * ldht]) {
                    result[11 - 1] = ulpinv;
                }
                if (vs[(i - 1) + (j - 1) * ldvs] != vs1[(i - 1) + (j - 1) * ldvs1]) {
                    result[12 - 1] = ulpinv;
                }
            }
        }
        if (sdim != sdim1) {
            result[13 - 1] = ulpinv;
        }
        //
        //        Compute RCONDE without VS, and compare
        //
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("N", &sort, _Rslect, "E", n, ht, lda, sdim1, wrt, wit, vs1, ldvs, rcnde1, rcndv1, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[14 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx6", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx6", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            goto statement_250;
        }
        //
        //        Perform test (14)
        //
        if (rcnde1 != rconde) {
            result[14 - 1] = ulpinv;
        }
        //
        //        Perform tests (10), (11), (12), and (13)
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                result[10 - 1] = ulpinv;
            }
            for (j = 1; j <= n; j = j + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != ht[(i - 1) + (j - 1) * ldht]) {
                    result[11 - 1] = ulpinv;
                }
                if (vs[(i - 1) + (j - 1) * ldvs] != vs1[(i - 1) + (j - 1) * ldvs1]) {
                    result[12 - 1] = ulpinv;
                }
            }
        }
        if (sdim != sdim1) {
            result[13 - 1] = ulpinv;
        }
        //
        //        Compute RCONDV with VS, and compare
        //
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("V", &sort, _Rslect, "V", n, ht, lda, sdim1, wrt, wit, vs1, ldvs, rcnde1, rcndv1, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[15 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx7", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx7", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            goto statement_250;
        }
        //
        //        Perform test (15)
        //
        if (rcndv1 != rcondv) {
            result[15 - 1] = ulpinv;
        }
        //
        //        Perform tests (10), (11), (12), and (13)
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                result[10 - 1] = ulpinv;
            }
            for (j = 1; j <= n; j = j + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != ht[(i - 1) + (j - 1) * ldht]) {
                    result[11 - 1] = ulpinv;
                }
                if (vs[(i - 1) + (j - 1) * ldvs] != vs1[(i - 1) + (j - 1) * ldvs1]) {
                    result[12 - 1] = ulpinv;
                }
            }
        }
        if (sdim != sdim1) {
            result[13 - 1] = ulpinv;
        }
        //
        //        Compute RCONDV without VS, and compare
        //
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("N", &sort, _Rslect, "V", n, ht, lda, sdim1, wrt, wit, vs1, ldvs, rcnde1, rcndv1, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[15 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeesx8", iinfo, n, jtype, iseed;
            } else {
                write(nounit, format_9999), "Rgeesx8", iinfo, n, iseed[1 - 1];
            }
            info = abs(iinfo);
            goto statement_250;
        }
        //
        //        Perform test (15)
        //
        if (rcndv1 != rcondv) {
            result[15 - 1] = ulpinv;
        }
        //
        //        Perform tests (10), (11), (12), and (13)
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                result[10 - 1] = ulpinv;
            }
            for (j = 1; j <= n; j = j + 1) {
                if (h[(i - 1) + (j - 1) * ldh] != ht[(i - 1) + (j - 1) * ldht]) {
                    result[11 - 1] = ulpinv;
                }
                if (vs[(i - 1) + (j - 1) * ldvs] != vs1[(i - 1) + (j - 1) * ldvs1]) {
                    result[12 - 1] = ulpinv;
                }
            }
        }
        if (sdim != sdim1) {
            result[13 - 1] = ulpinv;
        }
        //
    }
//
statement_250:
    //
    //     If there are precomputed reciprocal condition numbers, compare
    //     computed values with them.
    //
    if (comp) {
        //
        //        First set up SELOPT, SELDIM, SELVAL, SELWR, and SELWI so that
        //        the logical function _Rslect selects the eigenvalues specified
        //        by NSLCT and ISLCT.
        //
        seldim = n;
        selopt = 1;
        eps = max(ulp, epsin);
        for (i = 1; i <= n; i = i + 1) {
            ipnt[i - 1] = i;
            selval[i - 1] = false;
            selwr[i - 1] = wrtmp[i - 1];
            selwi[i - 1] = witmp[i - 1];
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            kmin = i;
            vrmin = wrtmp[i - 1];
            vimin = witmp[i - 1];
            for (j = i + 1; j <= n; j = j + 1) {
                if (wrtmp[j - 1] < vrmin) {
                    kmin = j;
                    vrmin = wrtmp[j - 1];
                    vimin = witmp[j - 1];
                }
            }
            wrtmp[kmin - 1] = wrtmp[i - 1];
            witmp[kmin - 1] = witmp[i - 1];
            wrtmp[i - 1] = vrmin;
            witmp[i - 1] = vimin;
            itmp = ipnt[i - 1];
            ipnt[i - 1] = ipnt[kmin - 1];
            ipnt[kmin - 1] = itmp;
        }
        for (i = 1; i <= nslct; i = i + 1) {
            selval[ipnt[islct[i - 1] - 1] - 1] = true;
        }
        //
        //        Compute condition numbers
        //
        Rlacpy("F", n, n, a, lda, ht, lda);
        Rgeesx("N", "S", _Rslect, "B", n, ht, lda, sdim1, wrt, wit, vs1, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, iinfo);
        if (iinfo != 0 && iinfo != n + 2) {
            result[16 - 1] = ulpinv;
            result[17 - 1] = ulpinv;
            write(nounit, format_9999), "Rgeesx9", iinfo, n, iseed[1 - 1];
            info = abs(iinfo);
            goto statement_300;
        }
        //
        //        Compare condition number for average of selected eigenvalues
        //        taking its condition number into account
        //
        anorm = Rlange("1", n, n, a, lda, work);
        v = max(castREAL(n) * eps * anorm, smlnum);
        if (anorm == zero) {
            v = one;
        }
        if (v > rcondv) {
            tol = one;
        } else {
            tol = v / rcondv;
        }
        if (v > rcdvin) {
            tolin = one;
        } else {
            tolin = v / rcdvin;
        }
        tol = max(tol, smlnum / eps);
        tolin = max(tolin, smlnum / eps);
        if (eps * (rcdein - tolin) > rconde + tol) {
            result[16 - 1] = ulpinv;
        } else if (rcdein - tolin > rconde + tol) {
            result[16 - 1] = (rcdein - tolin) / (rconde + tol);
        } else if (rcdein + tolin < eps * (rconde - tol)) {
            result[16 - 1] = ulpinv;
        } else if (rcdein + tolin < rconde - tol) {
            result[16 - 1] = (rconde - tol) / (rcdein + tolin);
        } else {
            result[16 - 1] = one;
        }
        //
        //        Compare condition numbers for right invariant subspace
        //        taking its condition number into account
        //
        if (v > rcondv * rconde) {
            tol = rcondv;
        } else {
            tol = v / rconde;
        }
        if (v > rcdvin * rcdein) {
            tolin = rcdvin;
        } else {
            tolin = v / rcdein;
        }
        tol = max(tol, smlnum / eps);
        tolin = max(tolin, smlnum / eps);
        if (eps * (rcdvin - tolin) > rcondv + tol) {
            result[17 - 1] = ulpinv;
        } else if (rcdvin - tolin > rcondv + tol) {
            result[17 - 1] = (rcdvin - tolin) / (rcondv + tol);
        } else if (rcdvin + tolin < eps * (rcondv - tol)) {
            result[17 - 1] = ulpinv;
        } else if (rcdvin + tolin < rcondv - tol) {
            result[17 - 1] = (rcondv - tol) / (rcdvin + tolin);
        } else {
            result[17 - 1] = one;
        }
    //
    statement_300:;
        //
    }
    //
    //     End of Rget24
    //
}

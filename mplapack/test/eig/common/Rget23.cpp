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

void Rget23(bool const comp, const char *balanc, INTEGER const jtype, REAL const thresh, INTEGER *iseed, INTEGER const nounit, INTEGER const n, REAL *a, INTEGER const lda, REAL *h, REAL *wr, REAL *wi, REAL *wr1, REAL *wi1, REAL *vl, INTEGER const ldvl, REAL *vr, INTEGER const ldvr, REAL *lre, INTEGER const ldlre, REAL *rcondv, REAL *rcndv1, REAL *rcdvin, REAL *rconde, REAL *rcnde1, REAL *rcdein, REAL *scale, REAL *scale1, REAL *result, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
    FEM_CMN_SVE(Rget23);
    iseed([4]);
    a([lda * star]);
    h([lda * star]);
    vl([ldvl * star]);
    vr([ldvr * star]);
    lre([ldlre * star]);
    result([11]);
    common_write write(cmn);
    str_arr_ref<1> sens(sve.sens, [2]);
    if (is_called_first_time) {
        static const char *values[] = {"N", "V"};
        data_of_type_str(FEM_VALUES_AND_SIZE), sens;
    }
    bool nobal = false;
    bool balok = false;
    const REAL zero = 0.0;
    INTEGER i = 0;
    const REAL one = 1.0;
    REAL ulp = 0.0;
    REAL smlnum = 0.0;
    REAL ulpinv = 0.0;
    char sense;
    INTEGER isensm = 0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL abnrm = 0.0;
    INTEGER iinfo = 0;
    REAL res[2];
    INTEGER j = 0;
    REAL tnrm = 0.0;
    REAL vmx = 0.0;
    REAL vrmx = 0.0;
    INTEGER jj = 0;
    REAL vtst = 0.0;
    const REAL two = 2.0;
    INTEGER isens = 0;
    REAL dum[1];
    INTEGER ilo1 = 0;
    INTEGER ihi1 = 0;
    REAL abnrm1 = 0.0;
    INTEGER kmin = 0;
    REAL vrmin = 0.0;
    REAL vimin = 0.0;
    const REAL epsin = 5.9605e-8;
    REAL eps = 0.0;
    REAL v = 0.0;
    REAL tol = 0.0;
    REAL tolin = 0.0;
    REAL vmax = 0.0;
    static const char *format_9998 = "(' Rget23: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                                     "', BALANC = ',a,', ISEED=(',3(i5,','),i5,')')";
    static const char *format_9999 = "(' Rget23: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
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
    //     Check for errors
    //
    nobal = Mlsame(balanc, "N");
    balok = nobal || Mlsame(balanc, "P") || Mlsame(balanc, "S") || Mlsame(balanc, "B");
    info = 0;
    if (!balok) {
        info = -2;
    } else if (thresh < zero) {
        info = -4;
    } else if (nounit <= 0) {
        info = -6;
    } else if (n < 0) {
        info = -7;
    } else if (lda < 1 || lda < n) {
        info = -9;
    } else if (ldvl < 1 || ldvl < n) {
        info = -16;
    } else if (ldvr < 1 || ldvr < n) {
        info = -18;
    } else if (ldlre < 1 || ldlre < n) {
        info = -20;
    } else if (lwork < 3 * n || (comp && lwork < 6 * n + n * n)) {
        info = -31;
    }
    //
    if (info != 0) {
        Mxerbla("Rget23", -info);
        return;
    }
    //
    //     Quick return if nothing to do
    //
    for (i = 1; i <= 11; i = i + 1) {
        result[i - 1] = -one;
    }
    //
    if (n == 0) {
        return;
    }
    //
    //     More Important constants
    //
    ulp = Rlamch("Precision");
    smlnum = Rlamch("S");
    ulpinv = one / ulp;
    //
    //     Compute eigenvalues and eigenvectors, and test them
    //
    if (lwork >= 6 * n + n * n) {
        sense = "B";
        isensm = 2;
    } else {
        sense = "E";
        isensm = 1;
    }
    Rlacpy("F", n, n, a, lda, h, lda);
    Rgeevx(balanc, "V", "V", sense, n, h, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, iinfo);
    if (iinfo != 0) {
        result[1 - 1] = ulpinv;
        if (jtype != 22) {
            write(nounit, format_9998), "Rgeevx1", iinfo, n, jtype, balanc, iseed;
        } else {
            write(nounit, format_9999), "Rgeevx1", iinfo, n, iseed(1);
        }
        info = abs(iinfo);
        return;
    }
    //
    //     Do Test (1)
    //
    Rget22("N", "N", "N", n, a, lda, vr, ldvr, wr, wi, work, res);
    result[1 - 1] = res[1 - 1];
    //
    //     Do Test (2)
    //
    Rget22("T", "N", "T", n, a, lda, vl, ldvl, wr, wi, work, res);
    result[2 - 1] = res[1 - 1];
    //
    //     Do Test (3)
    //
    for (j = 1; j <= n; j = j + 1) {
        tnrm = one;
        if (wi[j - 1] == zero) {
            tnrm = Rnrm2(n, &vr[(j - 1) * ldvr], 1);
        } else if (wi[j - 1] > zero) {
            tnrm = Rlapy2(Rnrm2(n, &vr[(j - 1) * ldvr], 1), Rnrm2(n, &vr[((j + 1) - 1) * ldvr], 1));
        }
        result[3 - 1] = max({result[3 - 1], min(ulpinv, abs(tnrm - one) / ulp)});
        if (wi[j - 1] > zero) {
            vmx = zero;
            vrmx = zero;
            for (jj = 1; jj <= n; jj = jj + 1) {
                vtst = Rlapy2(vr[(jj - 1) + (j - 1) * ldvr], &vr[(jj - 1) + ((j + 1) - 1) * ldvr]);
                if (vtst > vmx) {
                    vmx = vtst;
                }
                if (vr[(jj - 1) + ((j + 1) - 1) * ldvr] == zero && abs(vr[(jj - 1) + (j - 1) * ldvr]) > vrmx) {
                    vrmx = abs(vr[(jj - 1) + (j - 1) * ldvr]);
                }
            }
            if (vrmx / vmx < one - two * ulp) {
                result[3 - 1] = ulpinv;
            }
        }
    }
    //
    //     Do Test (4)
    //
    for (j = 1; j <= n; j = j + 1) {
        tnrm = one;
        if (wi[j - 1] == zero) {
            tnrm = Rnrm2(n, &vl[(j - 1) * ldvl], 1);
        } else if (wi[j - 1] > zero) {
            tnrm = Rlapy2(Rnrm2(n, &vl[(j - 1) * ldvl], 1), Rnrm2(n, &vl[((j + 1) - 1) * ldvl], 1));
        }
        result[4 - 1] = max({result[4 - 1], min(ulpinv, abs(tnrm - one) / ulp)});
        if (wi[j - 1] > zero) {
            vmx = zero;
            vrmx = zero;
            for (jj = 1; jj <= n; jj = jj + 1) {
                vtst = Rlapy2(vl[(jj - 1) + (j - 1) * ldvl], &vl[(jj - 1) + ((j + 1) - 1) * ldvl]);
                if (vtst > vmx) {
                    vmx = vtst;
                }
                if (vl[(jj - 1) + ((j + 1) - 1) * ldvl] == zero && abs(vl[(jj - 1) + (j - 1) * ldvl]) > vrmx) {
                    vrmx = abs(vl[(jj - 1) + (j - 1) * ldvl]);
                }
            }
            if (vrmx / vmx < one - two * ulp) {
                result[4 - 1] = ulpinv;
            }
        }
    }
    //
    //     Test for all options of computing condition numbers
    //
    for (isens = 1; isens <= isensm; isens = isens + 1) {
        //
        sense = sens[isens - 1];
        //
        //        Compute eigenvalues only, and test them
        //
        Rlacpy("F", n, n, a, lda, h, lda);
        Rgeevx(balanc, "N", "N", sense, n, h, lda, wr1, wi1, dum, 1, dum, 1, ilo1, ihi1, scale1, abnrm1, rcnde1, rcndv1, work, lwork, iwork, iinfo);
        if (iinfo != 0) {
            result[1 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeevx2", iinfo, n, jtype, balanc, iseed;
            } else {
                write(nounit, format_9999), "Rgeevx2", iinfo, n, iseed(1);
            }
            info = abs(iinfo);
            goto statement_190;
        }
        //
        //        Do Test (5)
        //
        for (j = 1; j <= n; j = j + 1) {
            if (wr[j - 1] != wr1[j - 1] || wi[j - 1] != wi1[j - 1]) {
                result[5 - 1] = ulpinv;
            }
        }
        //
        //        Do Test (8)
        //
        if (!nobal) {
            for (j = 1; j <= n; j = j + 1) {
                if (scale[j - 1] != scale1[j - 1]) {
                    result[8 - 1] = ulpinv;
                }
            }
            if (ilo != ilo1) {
                result[8 - 1] = ulpinv;
            }
            if (ihi != ihi1) {
                result[8 - 1] = ulpinv;
            }
            if (abnrm != abnrm1) {
                result[8 - 1] = ulpinv;
            }
        }
        //
        //        Do Test (9)
        //
        if (isens == 2 && n > 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (rcondv[j - 1] != rcndv1[j - 1]) {
                    result[9 - 1] = ulpinv;
                }
            }
        }
        //
        //        Compute eigenvalues and right eigenvectors, and test them
        //
        Rlacpy("F", n, n, a, lda, h, lda);
        Rgeevx(balanc, "N", "V", sense, n, h, lda, wr1, wi1, dum, 1, lre, ldlre, ilo1, ihi1, scale1, abnrm1, rcnde1, rcndv1, work, lwork, iwork, iinfo);
        if (iinfo != 0) {
            result[1 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeevx3", iinfo, n, jtype, balanc, iseed;
            } else {
                write(nounit, format_9999), "Rgeevx3", iinfo, n, iseed(1);
            }
            info = abs(iinfo);
            goto statement_190;
        }
        //
        //        Do Test (5) again
        //
        for (j = 1; j <= n; j = j + 1) {
            if (wr[j - 1] != wr1[j - 1] || wi[j - 1] != wi1[j - 1]) {
                result[5 - 1] = ulpinv;
            }
        }
        //
        //        Do Test (6)
        //
        for (j = 1; j <= n; j = j + 1) {
            for (jj = 1; jj <= n; jj = jj + 1) {
                if (vr[(j - 1) + (jj - 1) * ldvr] != lre[(j - 1) + (jj - 1) * ldlre]) {
                    result[6 - 1] = ulpinv;
                }
            }
        }
        //
        //        Do Test (8) again
        //
        if (!nobal) {
            for (j = 1; j <= n; j = j + 1) {
                if (scale[j - 1] != scale1[j - 1]) {
                    result[8 - 1] = ulpinv;
                }
            }
            if (ilo != ilo1) {
                result[8 - 1] = ulpinv;
            }
            if (ihi != ihi1) {
                result[8 - 1] = ulpinv;
            }
            if (abnrm != abnrm1) {
                result[8 - 1] = ulpinv;
            }
        }
        //
        //        Do Test (9) again
        //
        if (isens == 2 && n > 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (rcondv[j - 1] != rcndv1[j - 1]) {
                    result[9 - 1] = ulpinv;
                }
            }
        }
        //
        //        Compute eigenvalues and left eigenvectors, and test them
        //
        Rlacpy("F", n, n, a, lda, h, lda);
        Rgeevx(balanc, "V", "N", sense, n, h, lda, wr1, wi1, lre, ldlre, dum, 1, ilo1, ihi1, scale1, abnrm1, rcnde1, rcndv1, work, lwork, iwork, iinfo);
        if (iinfo != 0) {
            result[1 - 1] = ulpinv;
            if (jtype != 22) {
                write(nounit, format_9998), "Rgeevx4", iinfo, n, jtype, balanc, iseed;
            } else {
                write(nounit, format_9999), "Rgeevx4", iinfo, n, iseed(1);
            }
            info = abs(iinfo);
            goto statement_190;
        }
        //
        //        Do Test (5) again
        //
        for (j = 1; j <= n; j = j + 1) {
            if (wr[j - 1] != wr1[j - 1] || wi[j - 1] != wi1[j - 1]) {
                result[5 - 1] = ulpinv;
            }
        }
        //
        //        Do Test (7)
        //
        for (j = 1; j <= n; j = j + 1) {
            for (jj = 1; jj <= n; jj = jj + 1) {
                if (vl[(j - 1) + (jj - 1) * ldvl] != lre[(j - 1) + (jj - 1) * ldlre]) {
                    result[7 - 1] = ulpinv;
                }
            }
        }
        //
        //        Do Test (8) again
        //
        if (!nobal) {
            for (j = 1; j <= n; j = j + 1) {
                if (scale[j - 1] != scale1[j - 1]) {
                    result[8 - 1] = ulpinv;
                }
            }
            if (ilo != ilo1) {
                result[8 - 1] = ulpinv;
            }
            if (ihi != ihi1) {
                result[8 - 1] = ulpinv;
            }
            if (abnrm != abnrm1) {
                result[8 - 1] = ulpinv;
            }
        }
        //
        //        Do Test (9) again
        //
        if (isens == 2 && n > 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (rcondv[j - 1] != rcndv1[j - 1]) {
                    result[9 - 1] = ulpinv;
                }
            }
        }
    //
    statement_190:;
        //
    }
    //
    //     If COMP, compare condition numbers to precomputed ones
    //
    if (comp) {
        Rlacpy("F", n, n, a, lda, h, lda);
        Rgeevx("N", "V", "V", "B", n, h, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, iinfo);
        if (iinfo != 0) {
            result[1 - 1] = ulpinv;
            write(nounit, format_9999), "Rgeevx5", iinfo, n, iseed(1);
            info = abs(iinfo);
            goto statement_250;
        }
        //
        //        Sort eigenvalues and condition numbers lexicographically
        //        to compare with inputs
        //
        for (i = 1; i <= n - 1; i = i + 1) {
            kmin = i;
            vrmin = wr[i - 1];
            vimin = wi[i - 1];
            for (j = i + 1; j <= n; j = j + 1) {
                if (wr[j - 1] < vrmin) {
                    kmin = j;
                    vrmin = wr[j - 1];
                    vimin = wi[j - 1];
                }
            }
            wr[kmin - 1] = wr[i - 1];
            wi[kmin - 1] = wi[i - 1];
            wr[i - 1] = vrmin;
            wi[i - 1] = vimin;
            vrmin = rconde[kmin - 1];
            rconde[kmin - 1] = rconde[i - 1];
            rconde[i - 1] = vrmin;
            vrmin = rcondv[kmin - 1];
            rcondv[kmin - 1] = rcondv[i - 1];
            rcondv[i - 1] = vrmin;
        }
        //
        //        Compare condition numbers for eigenvectors
        //        taking their condition numbers into account
        //
        result[10 - 1] = zero;
        eps = max(epsin, ulp);
        v = max(n.real() * eps * abnrm, smlnum);
        if (abnrm == zero) {
            v = one;
        }
        for (i = 1; i <= n; i = i + 1) {
            if (v > rcondv[i - 1] * rconde[i - 1]) {
                tol = rcondv[i - 1];
            } else {
                tol = v / rconde[i - 1];
            }
            if (v > rcdvin[i - 1] * rcdein[i - 1]) {
                tolin = rcdvin[i - 1];
            } else {
                tolin = v / rcdein[i - 1];
            }
            tol = max(tol, smlnum / eps);
            tolin = max(tolin, smlnum / eps);
            if (eps * (rcdvin[i - 1] - tolin) > rcondv[i - 1] + tol) {
                vmax = one / eps;
            } else if (rcdvin[i - 1] - tolin > rcondv[i - 1] + tol) {
                vmax = (rcdvin[i - 1] - tolin) / (rcondv[i - 1] + tol);
            } else if (rcdvin[i - 1] + tolin < eps * (rcondv[i - 1] - tol)) {
                vmax = one / eps;
            } else if (rcdvin[i - 1] + tolin < rcondv[i - 1] - tol) {
                vmax = (rcondv[i - 1] - tol) / (rcdvin[i - 1] + tolin);
            } else {
                vmax = one;
            }
            result[10 - 1] = max(result[10 - 1], vmax);
        }
        //
        //        Compare condition numbers for eigenvalues
        //        taking their condition numbers into account
        //
        result[11 - 1] = zero;
        for (i = 1; i <= n; i = i + 1) {
            if (v > rcondv[i - 1]) {
                tol = one;
            } else {
                tol = v / rcondv[i - 1];
            }
            if (v > rcdvin[i - 1]) {
                tolin = one;
            } else {
                tolin = v / rcdvin[i - 1];
            }
            tol = max(tol, smlnum / eps);
            tolin = max(tolin, smlnum / eps);
            if (eps * (rcdein[i - 1] - tolin) > rconde[i - 1] + tol) {
                vmax = one / eps;
            } else if (rcdein[i - 1] - tolin > rconde[i - 1] + tol) {
                vmax = (rcdein[i - 1] - tolin) / (rconde[i - 1] + tol);
            } else if (rcdein[i - 1] + tolin < eps * (rconde[i - 1] - tol)) {
                vmax = one / eps;
            } else if (rcdein[i - 1] + tolin < rconde[i - 1] - tol) {
                vmax = (rconde[i - 1] - tol) / (rcdein[i - 1] + tolin);
            } else {
                vmax = one;
            }
            result[11 - 1] = max(result[11 - 1], vmax);
        }
    statement_250:;
        //
    }
    //
    //     End of Rget23
    //
}

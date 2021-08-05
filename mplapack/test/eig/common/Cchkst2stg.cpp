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

void Cchkst2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, COMPLEX *a, INTEGER const lda, COMPLEX *ap, REAL *sd, REAL *se, REAL *d1, REAL *d2, REAL *d3, REAL *d4, REAL *d5, REAL *wa1, REAL *wa2, REAL *wa3, REAL *wr, COMPLEX *u, INTEGER const ldu, COMPLEX *v, COMPLEX *vp, COMPLEX *tau, COMPLEX *z, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER *iwork, INTEGER const liwork, REAL *result, INTEGER &info) {
    INTEGER ldv = ldu;
    INTEGER ldz = ldu;
    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 21;
    INTEGER ktype[21] = {1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 10};
    INTEGER kmagn[21] = {1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 2, 3, 1};
    INTEGER kmode[21] = {0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 3, 1, 4, 4, 3};
    INTEGER idumma[1];
    INTEGER ntestt = 0;
    bool badnn = false;
    bool tryrac = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    INTEGER nblock = 0;
    REAL unfl = 0.0;
    const REAL one = 1.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    REAL ulpinv = 0.0;
    const REAL two = 2.0;
    INTEGER log2ui = 0;
    REAL rtunfl = 0.0;
    REAL rtovfl = 0.0;
    INTEGER i = 0;
    INTEGER iseed2[4];
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER lgn = 0;
    INTEGER lwedc = 0;
    INTEGER lrwedc = 0;
    INTEGER liwedc = 0;
    INTEGER nap = 0;
    REAL aninv = 0.0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ntest = 0;
    INTEGER ioldsd[4];
    INTEGER itype = 0;
    INTEGER imode = 0;
    REAL anorm = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    const REAL ten = 10.0;
    INTEGER jc = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const REAL zero = 0.0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    const REAL half = one / two;
    INTEGER lh = 0;
    INTEGER lw = 0;
    REAL temp3 = 0.0;
    REAL temp4 = 0.0;
    INTEGER jr = 0;
    REAL dumma[1];
    const REAL hun = 100.0;
    REAL vl = 0.0;
    REAL vu = 0.0;
    INTEGER il = 0;
    INTEGER iu = 0;
    REAL abstol = 0.0;
    INTEGER m = 0;
    INTEGER nsplit = 0;
    const REAL eight = 8.0;
    INTEGER itemp = 0;
    INTEGER m2 = 0;
    INTEGER m3 = 0;
    INTEGER inde = 0;
    INTEGER indrwk = 0;
    char buf[1024];
    const bool crel = false;
    const bool crange = false;
    static const char *format_9999 = "(' Cchkst2stg: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',3(i5,','),i5,')')";
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
    //     Keep ftnchek happy
    idumma[1 - 1] = 1;
    //
    //     Check for errors
    //
    ntestt = 0;
    info = 0;
    //
    //     Important constants
    //
    badnn = false;
    tryrac = true;
    nmax = 1;
    for (j = 1; j <= nsizes; j = j + 1) {
        nmax = max(nmax, nn[j - 1]);
        if (nn[j - 1] < 0) {
            badnn = true;
        }
    }
    //
    nblock = iMlaenv(1, "Chetrd", "L", nmax, -1, -1, -1);
    nblock = min({nmax, max((INTEGER)1, nblock)});
    //
    //     Check for errors
    //
    if (nsizes < 0) {
        info = -1;
    } else if (badnn) {
        info = -2;
    } else if (ntypes < 0) {
        info = -3;
    } else if (lda < nmax) {
        info = -9;
    } else if (ldu < nmax) {
        info = -23;
    } else if (2 * max((INTEGER)2, nmax) * max((INTEGER)2, nmax) > lwork) {
        info = -29;
    }
    //
    if (info != 0) {
        Mxerbla("Cchkst2stg", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (nsizes == 0 || ntypes == 0) {
        return;
    }
    //
    //     More Important constants
    //
    unfl = Rlamch("Safe minimum");
    ovfl = one / unfl;
    Rlabad(unfl, ovfl);
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    ulpinv = one / ulp;
    log2ui = castINTEGER(log(ulpinv) / log(two));
    rtunfl = sqrt(unfl);
    rtovfl = sqrt(ovfl);
    //
    //     Loop over sizes, types
    //
    for (i = 1; i <= 4; i = i + 1) {
        iseed2[i - 1] = iseed[i - 1];
    }
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        if (n > 0) {
            lgn = castINTEGER(log(castREAL(n)) / log(two));
            if ((INTEGER)pow((double)2, (double)lgn) < n) {
                lgn++;
            }
            if ((INTEGER)pow((double)2, (double)lgn) < n) {
                lgn++;
            }
            lwedc = 1 + 4 * n + 2 * n * lgn + 4 * n * n;
            lrwedc = 1 + 3 * n + 2 * n * lgn + 4 * n * n;
            liwedc = 6 + 6 * n + 5 * n * lgn;
        } else {
            lwedc = 8;
            lrwedc = 7;
            liwedc = 12;
        }
        nap = (n * (n + 1)) / 2;
        aninv = one / castREAL(max((INTEGER)1, n));
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_300;
            }
            nmats++;
            ntest = 0;
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            //           Compute "A"
            //
            //           Control parameters:
            //
            //               KMAGN  KMODE        KTYPE
            //           =1  O(1)   clustered 1  zero
            //           =2  large  clustered 2  identity
            //           =3  small  exponential  (none)
            //           =4         arithmetic   diagonal, (w/ eigenvalues)
            //           =5         random log   Hermitian, w/ eigenvalues
            //           =6         random       (none)
            //           =7                      random diagonal
            //           =8                      random Hermitian
            //           =9                      positive definite
            //           =10                     diagonally dominant tridiagonal
            //
            if (mtypes > maxtyp) {
                goto statement_100;
            }
            //
            itype = ktype[jtype - 1];
            imode = kmode[jtype - 1];
            //
            //           Compute norm
            //
            switch (kmagn[jtype - 1]) {
            case 1:
                goto statement_40;
            case 2:
                goto statement_50;
            case 3:
                goto statement_60;
            default:
                break;
            }
        //
        statement_40:
            anorm = one;
            goto statement_70;
        //
        statement_50:
            anorm = (rtovfl * ulp) * aninv;
            goto statement_70;
        //
        statement_60:
            anorm = rtunfl * n * ulpinv;
            goto statement_70;
        //
        statement_70:
            //
            Claset("Full", lda, n, czero, czero, a, lda);
            iinfo = 0;
            if (jtype <= 15) {
                cond = ulpinv;
            } else {
                cond = ulpinv * aninv / ten;
            }
            //
            //           Special Matrices -- Identity & Jordan block
            //
            //              Zero
            //
            if (itype == 1) {
                iinfo = 0;
                //
            } else if (itype == 2) {
                //
                //              Identity
                //
                for (jc = 1; jc <= n; jc = jc + 1) {
                    a[(jc - 1) + (jc - 1) * lda] = anorm;
                }
                //
            } else if (itype == 4) {
                //
                //              Diagonal Matrix, [Eigen]values Specified
                //
                Clatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, 0, 0, "N", a, lda, work, iinfo);
                //
            } else if (itype == 5) {
                //
                //              Hermitian, eigenvalues specified
                //
                Clatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, n, n, "N", a, lda, work, iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random eigenvalues
                //
                Clatmr(n, n, "S", iseed, "H", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Hermitian, random eigenvalues
                //
                Clatmr(n, n, "S", iseed, "H", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              Positive definite, eigenvalues specified.
                //
                Clatms(n, n, "S", iseed, "P", rwork, imode, cond, anorm, n, n, "N", a, lda, work, iinfo);
                //
            } else if (itype == 10) {
                //
                //              Positive definite tridiagonal, eigenvalues specified.
                //
                Clatms(n, n, "S", iseed, "P", rwork, imode, cond, anorm, 1, 1, "N", a, lda, work, iinfo);
                for (i = 2; i <= n; i = i + 1) {
                    temp1 = abs(a[((i - 1) - 1) + (i - 1) * lda]);
                    temp2 = sqrt(abs(a[((i - 1) - 1) + ((i - 1) - 1) * lda] * a[(i - 1) + (i - 1) * lda]));
                    if (temp1 > half * temp2) {
                        a[((i - 1) - 1) + (i - 1) * lda] = a[((i - 1) - 1) + (i - 1) * lda] * (half * temp2 / (unfl + temp1));
                        a[(i - 1) + ((i - 1) - 1) * lda] = conj(a[((i - 1) - 1) + (i - 1) * lda]);
                    }
                }
                //
            } else {
                //
                iinfo = 1;
            }
            //
            if (iinfo != 0) {
                write(nounit, format_9999), "Generator", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                return;
            }
        //
        statement_100:
            //
            //           Call Chetrd and Cungtr to compute S and U from
            //           upper triangle.
            //
            Clacpy("U", n, n, a, lda, v, ldu);
            //
            ntest = 1;
            Chetrd("U", n, v, ldu, sd, se, tau, work, lwork, iinfo);
            //
            if (iinfo != 0) {
                write(nounit, format_9999), "Chetrd(U)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[1 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            Clacpy("U", n, n, v, ldu, u, ldu);
            //
            ntest = 2;
            Cungtr("U", n, u, ldu, tau, work, lwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cungtr(U)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[2 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do tests 1 and 2
            //
            Chet21(2, "Upper", n, 1, a, lda, sd, se, u, ldu, v, ldu, tau, work, rwork, &result[1 - 1]);
            Chet21(3, "Upper", n, 1, a, lda, sd, se, u, ldu, v, ldu, tau, work, rwork, &result[2 - 1]);
            //
            //           Compute D1 the eigenvalues resulting from the tridiagonal
            //           form using the standard 1-stage algorithm and use it as a
            //           reference to compare with the 2-stage technique
            //
            //           Compute D1 from the 1-stage and used as reference for the
            //           2-stage
            //
            Rcopy(n, sd, 1, d1, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, rwork, 1);
            }
            //
            Csteqr("N", n, d1, rwork, work, ldu, &rwork[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Csteqr(N)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[3 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           2-STAGE TRD Upper case is used to compute D2.
            //           Note to set SD and SE to zero to be sure not reusing
            //           the one from above. Compare it with D1 computed
            //           using the 1-stage.
            //
            Rlaset("Full", n, 1, zero, zero, sd, n);
            Rlaset("Full", n, 1, zero, zero, se, n);
            Clacpy("U", n, n, a, lda, v, ldu);
            lh = max((INTEGER)1, 4 * n);
            lw = lwork - lh;
            Chetrd_2stage("N", "U", n, v, ldu, sd, se, tau, work, lh, &work[(lh + 1) - 1], lw, iinfo);
            //
            //           Compute D2 from the 2-stage Upper case
            //
            Rcopy(n, sd, 1, d2, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, rwork, 1);
            }
            //
            ntest = 3;
            Csteqr("N", n, d2, rwork, work, ldu, &rwork[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Csteqr(N)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[3 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           2-STAGE TRD Lower case is used to compute D3.
            //           Note to set SD and SE to zero to be sure not reusing
            //           the one from above. Compare it with D1 computed
            //           using the 1-stage.
            //
            Rlaset("Full", n, 1, zero, zero, sd, n);
            Rlaset("Full", n, 1, zero, zero, se, n);
            Clacpy("L", n, n, a, lda, v, ldu);
            Chetrd_2stage("N", "L", n, v, ldu, sd, se, tau, work, lh, &work[(lh + 1) - 1], lw, iinfo);
            //
            //           Compute D3 from the 2-stage Upper case
            //
            Rcopy(n, sd, 1, d3, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, rwork, 1);
            }
            //
            ntest = 4;
            Csteqr("N", n, d3, rwork, work, ldu, &rwork[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Csteqr(N)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[4 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do Tests 3 and 4 which are similar to 11 and 12 but with the
            //           D1 computed using the standard 1-stage reduction as reference
            //
            ntest = 4;
            temp1 = zero;
            temp2 = zero;
            temp3 = zero;
            temp4 = zero;
            //
            for (j = 1; j <= n; j = j + 1) {
                temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d2[j - 1]))});
                temp2 = max(temp2, REAL(abs(d1[j - 1] - d2[j - 1])));
                temp3 = max({temp3, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                temp4 = max(temp4, REAL(abs(d1[j - 1] - d3[j - 1])));
            }
            //
            result[3 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            result[4 - 1] = temp4 / max(unfl, REAL(ulp * max(temp3, temp4)));
            //
            //           Store the upper triangle of A in AP
            //
            i = 0;
            for (jc = 1; jc <= n; jc = jc + 1) {
                for (jr = 1; jr <= jc; jr = jr + 1) {
                    i++;
                    ap[i - 1] = a[(jr - 1) + (jc - 1) * lda];
                }
            }
            //
            //           Call Chptrd and Cupgtr to compute S and U from AP
            //
            Ccopy(nap, ap, 1, vp, 1);
            //
            ntest = 5;
            Chptrd("U", n, vp, sd, se, tau, iinfo);
            //
            if (iinfo != 0) {
                write(nounit, format_9999), "Chptrd(U)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[5 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            ntest = 6;
            Cupgtr("U", n, vp, tau, u, ldu, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cupgtr(U)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[6 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do tests 5 and 6
            //
            Chpt21(2, "Upper", n, 1, ap, sd, se, u, ldu, vp, tau, work, rwork, &result[5 - 1]);
            Chpt21(3, "Upper", n, 1, ap, sd, se, u, ldu, vp, tau, work, rwork, &result[6 - 1]);
            //
            //           Store the lower triangle of A in AP
            //
            i = 0;
            for (jc = 1; jc <= n; jc = jc + 1) {
                for (jr = jc; jr <= n; jr = jr + 1) {
                    i++;
                    ap[i - 1] = a[(jr - 1) + (jc - 1) * lda];
                }
            }
            //
            //           Call Chptrd and Cupgtr to compute S and U from AP
            //
            Ccopy(nap, ap, 1, vp, 1);
            //
            ntest = 7;
            Chptrd("L", n, vp, sd, se, tau, iinfo);
            //
            if (iinfo != 0) {
                write(nounit, format_9999), "Chptrd(L)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[7 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            ntest = 8;
            Cupgtr("L", n, vp, tau, u, ldu, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cupgtr(L)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[8 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            Chpt21(2, "Lower", n, 1, ap, sd, se, u, ldu, vp, tau, work, rwork, &result[7 - 1]);
            Chpt21(3, "Lower", n, 1, ap, sd, se, u, ldu, vp, tau, work, rwork, &result[8 - 1]);
            //
            //           Call Csteqr to compute D1, D2, and Z, do tests.
            //
            //           Compute D1 and Z
            //
            Rcopy(n, sd, 1, d1, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, rwork, 1);
            }
            Claset("Full", n, n, czero, cone, z, ldu);
            //
            ntest = 9;
            Csteqr("V", n, d1, rwork, z, ldu, &rwork[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Csteqr(V)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[9 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Compute D2
            //
            Rcopy(n, sd, 1, d2, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, rwork, 1);
            }
            //
            ntest = 11;
            Csteqr("N", n, d2, rwork, work, ldu, &rwork[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Csteqr(N)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[11 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Compute D3 (using PWK method)
            //
            Rcopy(n, sd, 1, d3, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, rwork, 1);
            }
            //
            ntest = 12;
            Rsterf(n, d3, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rsterf", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[12 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do Tests 9 and 10
            //
            Cstt21(n, 0, sd, se, d1, dumma, z, ldu, work, rwork, &result[9 - 1]);
            //
            //           Do Tests 11 and 12
            //
            temp1 = zero;
            temp2 = zero;
            temp3 = zero;
            temp4 = zero;
            //
            for (j = 1; j <= n; j = j + 1) {
                temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d2[j - 1]))});
                temp2 = max(temp2, REAL(abs(d1[j - 1] - d2[j - 1])));
                temp3 = max({temp3, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                temp4 = max(temp4, REAL(abs(d1[j - 1] - d3[j - 1])));
            }
            //
            result[11 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            result[12 - 1] = temp4 / max(unfl, REAL(ulp * max(temp3, temp4)));
            //
            //           Do Test 13 -- Sturm Sequence Test of Eigenvalues
            //                         Go up by factors of two until it succeeds
            //
            ntest = 13;
            temp1 = thresh * (half - ulp);
            //
            for (j = 0; j <= log2ui; j = j + 1) {
                Rstech(n, sd, se, d1, temp1, rwork, iinfo);
                if (iinfo == 0) {
                    goto statement_170;
                }
                temp1 = temp1 * two;
            }
        //
        statement_170:
            result[13 - 1] = temp1;
            //
            //           For positive definite matrices ( JTYPE.GT.15 ) call Cpteqr
            //           and do tests 14, 15, and 16 .
            //
            if (jtype > 15) {
                //
                //              Compute D4 and Z4
                //
                Rcopy(n, sd, 1, d4, 1);
                if (n > 0) {
                    Rcopy(n - 1, se, 1, rwork, 1);
                }
                Claset("Full", n, n, czero, cone, z, ldu);
                //
                ntest = 14;
                Cpteqr("V", n, d4, rwork, z, ldu, &rwork[(n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Cpteqr(V)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[14 - 1] = ulpinv;
                        goto statement_280;
                    }
                }
                //
                //              Do Tests 14 and 15
                //
                Cstt21(n, 0, sd, se, d4, dumma, z, ldu, work, rwork, &result[14 - 1]);
                //
                //              Compute D5
                //
                Rcopy(n, sd, 1, d5, 1);
                if (n > 0) {
                    Rcopy(n - 1, se, 1, rwork, 1);
                }
                //
                ntest = 16;
                Cpteqr("N", n, d5, rwork, z, ldu, &rwork[(n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Cpteqr(N)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[16 - 1] = ulpinv;
                        goto statement_280;
                    }
                }
                //
                //              Do Test 16
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d4[j - 1])), REAL(abs(d5[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d4[j - 1] - d5[j - 1])));
                }
                //
                result[16 - 1] = temp2 / max(unfl, REAL(hun * ulp * max(temp1, temp2)));
            } else {
                result[14 - 1] = zero;
                result[15 - 1] = zero;
                result[16 - 1] = zero;
            }
            //
            //           Call Rstebz with different options and do tests 17-18.
            //
            //              If S is positive definite and diagonally dominant,
            //              ask for all eigenvalues with high relative accuracy.
            //
            vl = zero;
            vu = zero;
            il = 0;
            iu = 0;
            if (jtype == 21) {
                ntest = 17;
                abstol = unfl + unfl;
                Rstebz("A", "E", n, vl, vu, il, iu, abstol, sd, se, m, nsplit, wr, &iwork[1 - 1], &iwork[(n + 1) - 1], rwork, &iwork[(2 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstebz(A,rel)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[17 - 1] = ulpinv;
                        goto statement_280;
                    }
                }
                //
                //              Do test 17
                //
                temp2 = two * (two * n - one) * ulp * (one + eight * half * half) / ((one - half) * (one - half) * (one - half) * (one - half));
                //
                temp1 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max(temp1, REAL(abs(d4[j - 1] - wr[(n - j + 1) - 1]) / (abstol + abs(d4[j - 1]))));
                }
                //
                result[17 - 1] = temp1 / temp2;
            } else {
                result[17 - 1] = zero;
            }
            //
            //           Now ask for all eigenvalues with high absolute accuracy.
            //
            ntest = 18;
            abstol = unfl + unfl;
            Rstebz("A", "E", n, vl, vu, il, iu, abstol, sd, se, m, nsplit, wa1, &iwork[1 - 1], &iwork[(n + 1) - 1], rwork, &iwork[(2 * n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rstebz(A)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[18 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do test 18
            //
            temp1 = zero;
            temp2 = zero;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = max({temp1, REAL(abs(d3[j - 1])), REAL(abs(wa1[j - 1]))});
                temp2 = max(temp2, REAL(abs(d3[j - 1] - wa1[j - 1])));
            }
            //
            result[18 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            //           Choose random values for IL and IU, and ask for the
            //           IL-th through IU-th eigenvalues.
            //
            ntest = 19;
            if (n <= 1) {
                il = 1;
                iu = n;
            } else {
                il = 1 + castINTEGER(castREAL(n - 1) * Rlarnd(1, iseed2));
                iu = 1 + castINTEGER(castREAL(n - 1) * Rlarnd(1, iseed2));
                if (iu < il) {
                    itemp = iu;
                    iu = il;
                    il = itemp;
                }
            }
            //
            Rstebz("I", "E", n, vl, vu, il, iu, abstol, sd, se, m2, nsplit, wa2, &iwork[1 - 1], &iwork[(n + 1) - 1], rwork, &iwork[(2 * n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rstebz(I)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[19 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Determine the values VL and VU of the IL-th and IU-th
            //           eigenvalues and ask for all eigenvalues in this range.
            //
            if (n > 0) {
                if (il != 1) {
                    vl = wa1[il - 1] - max({REAL(half * (wa1[il - 1] - wa1[(il - 1) - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                } else {
                    vl = wa1[1 - 1] - max({REAL(half * (wa1[n - 1] - wa1[1 - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                }
                if (iu != n) {
                    vu = wa1[iu - 1] + max({REAL(half * (wa1[(iu + 1) - 1] - wa1[iu - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                } else {
                    vu = wa1[n - 1] + max({REAL(half * (wa1[n - 1] - wa1[1 - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                }
            } else {
                vl = zero;
                vu = one;
            }
            //
            Rstebz("V", "E", n, vl, vu, il, iu, abstol, sd, se, m3, nsplit, wa3, &iwork[1 - 1], &iwork[(n + 1) - 1], rwork, &iwork[(2 * n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rstebz(V)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[19 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            if (m3 == 0 && n != 0) {
                result[19 - 1] = ulpinv;
                goto statement_280;
            }
            //
            //           Do test 19
            //
            temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
            temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
            if (n > 0) {
                temp3 = max(abs(wa1[n - 1]), abs(wa1[1 - 1]));
            } else {
                temp3 = zero;
            }
            //
            result[19 - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
            //
            //           Call Cstein to compute eigenvectors corresponding to
            //           eigenvalues in WA1.  (First call Rstebz again, to make sure
            //           it returns these eigenvalues in the correct order.)
            //
            ntest = 21;
            Rstebz("A", "B", n, vl, vu, il, iu, abstol, sd, se, m, nsplit, wa1, &iwork[1 - 1], &iwork[(n + 1) - 1], rwork, &iwork[(2 * n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rstebz(A,B)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[20 - 1] = ulpinv;
                    result[21 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            Cstein(n, sd, se, m, wa1, &iwork[1 - 1], &iwork[(n + 1) - 1], z, ldu, rwork, &iwork[(2 * n + 1) - 1], &iwork[(3 * n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cstein", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[20 - 1] = ulpinv;
                    result[21 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do tests 20 and 21
            //
            Cstt21(n, 0, sd, se, wa1, dumma, z, ldu, work, rwork, &result[20 - 1]);
            //
            //           Call Cstedc(I) to compute D1 and Z, do tests.
            //
            //           Compute D1 and Z
            //
            inde = 1;
            indrwk = inde + n;
            Rcopy(n, sd, 1, d1, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, &rwork[inde - 1], 1);
            }
            Claset("Full", n, n, czero, cone, z, ldu);
            //
            ntest = 22;
            Cstedc("I", n, d1, &rwork[inde - 1], z, ldu, work, lwedc, &rwork[indrwk - 1], lrwedc, iwork, liwedc, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cstedc(I)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[22 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do Tests 22 and 23
            //
            Cstt21(n, 0, sd, se, d1, dumma, z, ldu, work, rwork, &result[22 - 1]);
            //
            //           Call Cstedc(V) to compute D1 and Z, do tests.
            //
            //           Compute D1 and Z
            //
            Rcopy(n, sd, 1, d1, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, &rwork[inde - 1], 1);
            }
            Claset("Full", n, n, czero, cone, z, ldu);
            //
            ntest = 24;
            Cstedc("V", n, d1, &rwork[inde - 1], z, ldu, work, lwedc, &rwork[indrwk - 1], lrwedc, iwork, liwedc, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cstedc(V)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[24 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do Tests 24 and 25
            //
            Cstt21(n, 0, sd, se, d1, dumma, z, ldu, work, rwork, &result[24 - 1]);
            //
            //           Call Cstedc(N) to compute D2, do tests.
            //
            //           Compute D2
            //
            Rcopy(n, sd, 1, d2, 1);
            if (n > 0) {
                Rcopy(n - 1, se, 1, &rwork[inde - 1], 1);
            }
            Claset("Full", n, n, czero, cone, z, ldu);
            //
            ntest = 26;
            Cstedc("N", n, d2, &rwork[inde - 1], z, ldu, work, lwedc, &rwork[indrwk - 1], lrwedc, iwork, liwedc, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cstedc(N)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[26 - 1] = ulpinv;
                    goto statement_280;
                }
            }
            //
            //           Do Test 26
            //
            temp1 = zero;
            temp2 = zero;
            //
            for (j = 1; j <= n; j = j + 1) {
                temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d2[j - 1]))});
                temp2 = max(temp2, REAL(abs(d1[j - 1] - d2[j - 1])));
            }
            //
            result[26 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            //           Only test Cstemr if IEEE compliant
            //
            if (iMlaenv(10, "Cstemr", "VA", 1, 0, 0, 0) == 1 && iMlaenv(11, "Cstemr", "VA", 1, 0, 0, 0) == 1) {
                //
                //           Call Cstemr, do test 27 (relative eigenvalue accuracy)
                //
                //              If S is positive definite and diagonally dominant,
                //              ask for all eigenvalues with high relative accuracy.
                //
                vl = zero;
                vu = zero;
                il = 0;
                iu = 0;
                if (jtype == 21 && crel) {
                    ntest = 27;
                    abstol = unfl + unfl;
                    Cstemr("V", "A", n, sd, se, vl, vu, il, iu, m, wr, z, ldu, n, &iwork[1 - 1], tryrac, rwork, lrwork, &iwork[(2 * n + 1) - 1], lwork - 2 * n, iinfo);
                    if (iinfo != 0) {
                        write(nounit, format_9999), "Cstemr(V,A,rel)", iinfo, n, jtype, ioldsd;
                        info = abs(iinfo);
                        if (iinfo < 0) {
                            return;
                        } else {
                            result[27 - 1] = ulpinv;
                            goto statement_270;
                        }
                    }
                    //
                    //              Do test 27
                    //
                    temp2 = two * (two * n - one) * ulp * (one + eight * half * half) / ((one - half) * (one - half) * (one - half) * (one - half));
                    //
                    temp1 = zero;
                    for (j = 1; j <= n; j = j + 1) {
                        temp1 = max(temp1, REAL(abs(d4[j - 1] - wr[(n - j + 1) - 1]) / (abstol + abs(d4[j - 1]))));
                    }
                    //
                    result[27 - 1] = temp1 / temp2;
                    //
                    il = 1 + castINTEGER(castREAL(n - 1) * Rlarnd(1, iseed2));
                    iu = 1 + castINTEGER(castREAL(n - 1) * Rlarnd(1, iseed2));
                    if (iu < il) {
                        itemp = iu;
                        iu = il;
                        il = itemp;
                    }
                    //
                    if (crange) {
                        ntest = 28;
                        abstol = unfl + unfl;
                        Cstemr("V", "I", n, sd, se, vl, vu, il, iu, m, wr, z, ldu, n, &iwork[1 - 1], tryrac, rwork, lrwork, &iwork[(2 * n + 1) - 1], lwork - 2 * n, iinfo);
                        //
                        if (iinfo != 0) {
                            write(nounit, format_9999), "Cstemr(V,I,rel)", iinfo, n, jtype, ioldsd;
                            info = abs(iinfo);
                            if (iinfo < 0) {
                                return;
                            } else {
                                result[28 - 1] = ulpinv;
                                goto statement_270;
                            }
                        }
                        //
                        //                 Do test 28
                        //
                        temp2 = two * (two * n - one) * ulp * (one + eight * half * half) / ((one - half) * (one - half) * (one - half) * (one - half));
                        //
                        temp1 = zero;
                        for (j = il; j <= iu; j = j + 1) {
                            temp1 = max(temp1, REAL(abs(wr[(j - il + 1) - 1] - d4[(n - j + 1) - 1]) / (abstol + abs(wr[(j - il + 1) - 1]))));
                        }
                        //
                        result[28 - 1] = temp1 / temp2;
                    } else {
                        result[28 - 1] = zero;
                    }
                } else {
                    result[27 - 1] = zero;
                    result[28 - 1] = zero;
                }
                //
                //           Call Cstemr(V,I) to compute D1 and Z, do tests.
                //
                //           Compute D1 and Z
                //
                Rcopy(n, sd, 1, d5, 1);
                if (n > 0) {
                    Rcopy(n - 1, se, 1, rwork, 1);
                }
                Claset("Full", n, n, czero, cone, z, ldu);
                //
                if (crange) {
                    ntest = 29;
                    il = 1 + castINTEGER(castREAL(n - 1) * Rlarnd(1, iseed2));
                    iu = 1 + castINTEGER(castREAL(n - 1) * Rlarnd(1, iseed2));
                    if (iu < il) {
                        itemp = iu;
                        iu = il;
                        il = itemp;
                    }
                    Cstemr("V", "I", n, d5, rwork, vl, vu, il, iu, m, d1, z, ldu, n, &iwork[1 - 1], tryrac, &rwork[(n + 1) - 1], lrwork - n, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                    if (iinfo != 0) {
                        write(nounit, format_9999), "Cstemr(V,I)", iinfo, n, jtype, ioldsd;
                        info = abs(iinfo);
                        if (iinfo < 0) {
                            return;
                        } else {
                            result[29 - 1] = ulpinv;
                            goto statement_280;
                        }
                    }
                    //
                    //           Do Tests 29 and 30
                    //
                    //           Call Cstemr to compute D2, do tests.
                    //
                    //           Compute D2
                    //
                    Rcopy(n, sd, 1, d5, 1);
                    if (n > 0) {
                        Rcopy(n - 1, se, 1, rwork, 1);
                    }
                    //
                    ntest = 31;
                    Cstemr("N", "I", n, d5, rwork, vl, vu, il, iu, m, d2, z, ldu, n, &iwork[1 - 1], tryrac, &rwork[(n + 1) - 1], lrwork - n, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                    if (iinfo != 0) {
                        write(nounit, format_9999), "Cstemr(N,I)", iinfo, n, jtype, ioldsd;
                        info = abs(iinfo);
                        if (iinfo < 0) {
                            return;
                        } else {
                            result[31 - 1] = ulpinv;
                            goto statement_280;
                        }
                    }
                    //
                    //           Do Test 31
                    //
                    temp1 = zero;
                    temp2 = zero;
                    //
                    for (j = 1; j <= iu - il + 1; j = j + 1) {
                        temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d2[j - 1]))});
                        temp2 = max(temp2, REAL(abs(d1[j - 1] - d2[j - 1])));
                    }
                    //
                    result[31 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
                    //
                    //           Call Cstemr(V,V) to compute D1 and Z, do tests.
                    //
                    //           Compute D1 and Z
                    //
                    Rcopy(n, sd, 1, d5, 1);
                    if (n > 0) {
                        Rcopy(n - 1, se, 1, rwork, 1);
                    }
                    Claset("Full", n, n, czero, cone, z, ldu);
                    //
                    ntest = 32;
                    //
                    if (n > 0) {
                        if (il != 1) {
                            vl = d2[il - 1] - max({REAL(half * (d2[il - 1] - d2[(il - 1) - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                        } else {
                            vl = d2[1 - 1] - max({REAL(half * (d2[n - 1] - d2[1 - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                        }
                        if (iu != n) {
                            vu = d2[iu - 1] + max({REAL(half * (d2[(iu + 1) - 1] - d2[iu - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                        } else {
                            vu = d2[n - 1] + max({REAL(half * (d2[n - 1] - d2[1 - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                        }
                    } else {
                        vl = zero;
                        vu = one;
                    }
                    //
                    Cstemr("V", "V", n, d5, rwork, vl, vu, il, iu, m, d1, z, ldu, m, &iwork[1 - 1], tryrac, &rwork[(n + 1) - 1], lrwork - n, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                    if (iinfo != 0) {
                        write(nounit, format_9999), "Cstemr(V,V)", iinfo, n, jtype, ioldsd;
                        info = abs(iinfo);
                        if (iinfo < 0) {
                            return;
                        } else {
                            result[32 - 1] = ulpinv;
                            goto statement_280;
                        }
                    }
                    //
                    //           Do Tests 32 and 33
                    //
                    Cstt22(n, m, 0, sd, se, d1, dumma, z, ldu, work, m, rwork, &result[32 - 1]);
                    //
                    //           Call Cstemr to compute D2, do tests.
                    //
                    //           Compute D2
                    //
                    Rcopy(n, sd, 1, d5, 1);
                    if (n > 0) {
                        Rcopy(n - 1, se, 1, rwork, 1);
                    }
                    //
                    ntest = 34;
                    Cstemr("N", "V", n, d5, rwork, vl, vu, il, iu, m, d2, z, ldu, n, &iwork[1 - 1], tryrac, &rwork[(n + 1) - 1], lrwork - n, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                    if (iinfo != 0) {
                        write(nounit, format_9999), "Cstemr(N,V)", iinfo, n, jtype, ioldsd;
                        info = abs(iinfo);
                        if (iinfo < 0) {
                            return;
                        } else {
                            result[34 - 1] = ulpinv;
                            goto statement_280;
                        }
                    }
                    //
                    //           Do Test 34
                    //
                    temp1 = zero;
                    temp2 = zero;
                    //
                    for (j = 1; j <= iu - il + 1; j = j + 1) {
                        temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d2[j - 1]))});
                        temp2 = max(temp2, REAL(abs(d1[j - 1] - d2[j - 1])));
                    }
                    //
                    result[34 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
                } else {
                    result[29 - 1] = zero;
                    result[30 - 1] = zero;
                    result[31 - 1] = zero;
                    result[32 - 1] = zero;
                    result[33 - 1] = zero;
                    result[34 - 1] = zero;
                }
                //
                //           Call Cstemr(V,A) to compute D1 and Z, do tests.
                //
                //           Compute D1 and Z
                //
                Rcopy(n, sd, 1, d5, 1);
                if (n > 0) {
                    Rcopy(n - 1, se, 1, rwork, 1);
                }
                //
                ntest = 35;
                //
                Cstemr("V", "A", n, d5, rwork, vl, vu, il, iu, m, d1, z, ldu, n, &iwork[1 - 1], tryrac, &rwork[(n + 1) - 1], lrwork - n, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Cstemr(V,A)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[35 - 1] = ulpinv;
                        goto statement_280;
                    }
                }
                //
                //           Do Tests 35 and 36
                //
                Cstt22(n, m, 0, sd, se, d1, dumma, z, ldu, work, m, rwork, &result[35 - 1]);
                //
                //           Call Cstemr to compute D2, do tests.
                //
                //           Compute D2
                //
                Rcopy(n, sd, 1, d5, 1);
                if (n > 0) {
                    Rcopy(n - 1, se, 1, rwork, 1);
                }
                //
                ntest = 37;
                Cstemr("N", "A", n, d5, rwork, vl, vu, il, iu, m, d2, z, ldu, n, &iwork[1 - 1], tryrac, &rwork[(n + 1) - 1], lrwork - n, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Cstemr(N,A)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[37 - 1] = ulpinv;
                        goto statement_280;
                    }
                }
                //
                //           Do Test 34
                //
                temp1 = zero;
                temp2 = zero;
                //
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d2[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d1[j - 1] - d2[j - 1])));
                }
                //
                result[37 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            }
        statement_270:
        statement_280:
            ntestt += ntest;
            //
            //           End of Loop -- Check for RESULT(j) > THRESH
            //
            //           Print out tests which fail.
            //
            for (jr = 1; jr <= ntest; jr = jr + 1) {
                if (result[jr - 1] >= thresh) {
                    //
                    //                 If this is the first test to fail,
                    //                 print a header to the data file.
                    //
                    if (nerrs == 0) {
                        write(nounit, "(/,1x,a3,' -- Complex Hermitian eigenvalue problem')"), "ZST";
                        write(nounit, "(' Matrix types (see Cchkst2stg for details): ')");
                        write(nounit, "(/,' Special Matrices:',/,"
                                      "'  1=Zero matrix.                        ',"
                                      "'  5=Diagonal: clustered entries.',/,"
                                      "'  2=Identity matrix.                    ',"
                                      "'  6=Diagonal: large, evenly spaced.',/,"
                                      "'  3=Diagonal: evenly spaced entries.    ',"
                                      "'  7=Diagonal: small, evenly spaced.',/,"
                                      "'  4=Diagonal: geometr. spaced entries.')");
                        write(nounit, "(' Dense ',a,' Matrices:',/,"
                                      "'  8=Evenly spaced eigenvals.            ',"
                                      "' 12=Small, evenly spaced eigenvals.',/,"
                                      "'  9=Geometrically spaced eigenvals.     ',"
                                      "' 13=Matrix with random O(1) entries.',/,"
                                      "' 10=Clustered eigenvalues.              ',"
                                      "' 14=Matrix with large random entries.',/,"
                                      "' 11=Large, evenly spaced eigenvals.     ',"
                                      "' 15=Matrix with small random entries.')"),
                            "Hermitian";
                        write(nounit, "(' 16=Positive definite, evenly spaced eigenvalues',/,"
                                      "' 17=Positive definite, geometrically spaced eigenvlaues',/,"
                                      "' 18=Positive definite, clustered eigenvalues',/,"
                                      "' 19=Positive definite, small evenly spaced eigenvalues',/,"
                                      "' 20=Positive definite, large evenly spaced eigenvalues',/,"
                                      "' 21=Diagonally dominant tridiagonal, geometrically',"
                                      "' spaced eigenvalues')");
                        //
                        //                    Tests performed
                        //
                        write(nounit, "(/,'Test performed:  see Cchkst2stg for details.',/)");
                    }
                    nerrs++;
                    if (result[jr - 1] < 10000.0) {
                        sprintnum_short(buf, result[jr - 1]);
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i3,' is',0p,a)"),
                            n, jtype, ioldsd, jr, buf;
                    } else {
                        sprintnum_short(buf, result[jr - 1]);
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i3,' is',1p,a)"),
                            n, jtype, ioldsd, jr, buf;
                    }
                }
            }
        statement_300:;
        }
    }
    //
    //     Summary
    //
    Rlasum("ZST", nounit, nerrs, ntestt);
    //
    //     End of Cchkst2stg
    //
}

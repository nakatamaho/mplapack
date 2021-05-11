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

void Rdrvst2stg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *d1, REAL *d2, REAL *d3, REAL *d4, REAL *eveigs, REAL *wa1, REAL *wa2, REAL *wa3, REAL *u, INTEGER const ldu, REAL *v, REAL *tau, REAL *z, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, REAL *result, INTEGER &info) {
    INTEGER ldv = ldu;
    INTEGER ldz = ldu;
    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 18;
    INTEGER ktype[18] = {1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9};
    INTEGER kmagn[18] = {1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3};
    INTEGER kmode[18] = {0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4};
    const REAL zero = 0.0;
    REAL vl = 0.0;
    REAL vu = 0.0;
    INTEGER ntestt = 0;
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    REAL unfl = 0.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL ulpinv = 0.0;
    REAL rtunfl = 0.0;
    REAL rtovfl = 0.0;
    INTEGER i = 0;
    INTEGER iseed2[4];
    INTEGER iseed3[4];
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    const REAL two = 2.0;
    INTEGER lgn = 0;
    INTEGER lwedc = 0;
    INTEGER liwedc = 0;
    REAL aninv = 0.0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ntest = 0;
    INTEGER ioldsd[4];
    INTEGER itype = 0;
    INTEGER imode = 0;
    REAL anorm = 0.0;
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    INTEGER jcol = 0;
    INTEGER idumma[1];
    INTEGER ihbw = 0;
    INTEGER idiag = 0;
    INTEGER irow = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    REAL abstol = 0.0;
    INTEGER il = 0;
    INTEGER iu = 0;
    INTEGER itemp = 0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    INTEGER m = 0;
    REAL temp3 = 0.0;
    INTEGER m2 = 0;
    INTEGER m3 = 0;
    const REAL half = 0.5e0;
    const REAL ten = 10.0;
    INTEGER iuplo = 0;
    char uplo;
    INTEGER indx = 0;
    INTEGER kd = 0;
    static const char *format_9999 = "(' Rdrvst2stg: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Keep ftrnchek happy
    //
    vl = zero;
    vu = zero;
    //
    //     1)      Check for errors
    //
    ntestt = 0;
    info = 0;
    //
    badnn = false;
    nmax = 1;
    for (j = 1; j <= nsizes; j = j + 1) {
        nmax = max(nmax, nn[j - 1]);
        if (nn[j - 1] < 0) {
            badnn = true;
        }
    }
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
        info = -16;
    } else if (2 * pow2(max((INTEGER)2, nmax)) > lwork) {
        info = -21;
    }
    //
    if (info != 0) {
        Mxerbla("Rdrvst2stg", -info);
        return;
    }
    //
    //     Quick return if nothing to do
    //
    if (nsizes == 0 || ntypes == 0) {
        return;
    }
    //
    //     More Important constants
    //
    unfl = Rlamch("Safe minimum");
    ovfl = Rlamch("Overflow");
    Rlabad(unfl, ovfl);
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    ulpinv = one / ulp;
    rtunfl = sqrt(unfl);
    rtovfl = sqrt(ovfl);
    //
    //     Loop over sizes, types
    //
    for (i = 1; i <= 4; i = i + 1) {
        iseed2[i - 1] = iseed[i - 1];
        iseed3[i - 1] = iseed[i - 1];
    }
    //
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        if (n > 0) {
            lgn = castINTEGER(log(castREAL(n)) / log(two));
            if (pow(2, lgn) < n) {
                lgn++;
            }
            if (pow(2, lgn) < n) {
                lgn++;
            }
            lwedc = 1 + 4 * n + 2 * n * lgn + 4 * pow2(n);
            //           LIWEDC = 6 + 6*N + 5*N*LGN
            liwedc = 3 + 5 * n;
        } else {
            lwedc = 9;
            //           LIWEDC = 12
            liwedc = 8;
        }
        aninv = one / castREAL(max((INTEGER)1, n));
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            //
            if (!dotype[jtype - 1]) {
                goto statement_1730;
            }
            nmats++;
            ntest = 0;
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            //           2)      Compute "A"
            //
            //                   Control parameters:
            //
            //               KMAGN  KMODE        KTYPE
            //           =1  O(1)   clustered 1  zero
            //           =2  large  clustered 2  identity
            //           =3  small  exponential  (none)
            //           =4         arithmetic   diagonal, (w/ eigenvalues)
            //           =5         random log   symmetric, w/ eigenvalues
            //           =6         random       (none)
            //           =7                      random diagonal
            //           =8                      random symmetric
            //           =9                      band symmetric, w/ eigenvalues
            //
            if (mtypes > maxtyp) {
                goto statement_110;
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
            Rlaset("Full", lda, n, zero, zero, a, lda);
            iinfo = 0;
            cond = ulpinv;
            //
            //           Special Matrices -- Identity & Jordan block
            //
            //                   Zero
            //
            if (itype == 1) {
                iinfo = 0;
                //
            } else if (itype == 2) {
                //
                //              Identity
                //
                for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                    a[(jcol - 1) + (jcol - 1) * lda] = anorm;
                }
                //
            } else if (itype == 4) {
                //
                //              Diagonal Matrix, [Eigen]values Specified
                //
                Rlatms(n, n, "S", iseed, "S", work, imode, cond, anorm, 0, 0, "N", a, lda, &work[(n + 1) - 1], iinfo);
                //
            } else if (itype == 5) {
                //
                //              Symmetric, eigenvalues specified
                //
                Rlatms(n, n, "S", iseed, "S", work, imode, cond, anorm, n, n, "N", a, lda, &work[(n + 1) - 1], iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random eigenvalues
                //
                idumma[1 - 1] = 1;
                Rlatmr(n, n, "S", iseed, "S", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Symmetric, random eigenvalues
                //
                idumma[1 - 1] = 1;
                Rlatmr(n, n, "S", iseed, "S", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              Symmetric banded, eigenvalues specified
                //
                ihbw = castINTEGER((n - 1) * Rlarnd(1, iseed3));
                Rlatms(n, n, "S", iseed, "S", work, imode, cond, anorm, ihbw, ihbw, "Z", u, ldu, &work[(n + 1) - 1], iinfo);
                //
                //              Store as dense matrix for most routines.
                //
                Rlaset("Full", lda, n, zero, zero, a, lda);
                for (idiag = -ihbw; idiag <= ihbw; idiag = idiag + 1) {
                    irow = ihbw - idiag + 1;
                    j1 = max((INTEGER)1, idiag + 1);
                    j2 = min(n, n + idiag);
                    for (j = j1; j <= j2; j = j + 1) {
                        i = j - idiag;
                        a[(i - 1) + (j - 1) * lda] = u[(irow - 1) + (j - 1) * ldu];
                    }
                }
            } else {
                iinfo = 1;
            }
            //
            if (iinfo != 0) {
                write(nounit, format_9999), "Generator", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                return;
            }
        //
        statement_110:
            //
            abstol = unfl + unfl;
            if (n <= 1) {
                il = 1;
                iu = n;
            } else {
                il = 1 + castINTEGER((n - 1) * Rlarnd(1, iseed2));
                iu = 1 + castINTEGER((n - 1) * Rlarnd(1, iseed2));
                if (il > iu) {
                    itemp = il;
                    il = iu;
                    iu = itemp;
                }
            }
            //
            //           3)      If matrix is tridiagonal, call Rstev and Rstevx.
            //
            if (jtype <= 7) {
                ntest = 1;
                for (i = 1; i <= n; i = i + 1) {
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstev("V", n, d1, d2, z, ldu, work, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstev(V)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[1 - 1] = ulpinv;
                        result[2 - 1] = ulpinv;
                        result[3 - 1] = ulpinv;
                        goto statement_180;
                    }
                }
                //
                //              Do tests 1 and 2.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt21(n, 0, d3, d4, d1, d2, z, ldu, work, &result[1 - 1]);
                //
                ntest = 3;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstev("N", n, d3, d4, z, ldu, work, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstev(N)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[3 - 1] = ulpinv;
                        goto statement_180;
                    }
                }
                //
                //              Do test 3.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(d1[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(d1[j - 1] - d3[j - 1]));
                }
                result[3 - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_180:
                //
                ntest = 4;
                for (i = 1; i <= n; i = i + 1) {
                    eveigs[i - 1] = d3[i - 1];
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevx("V", "A", n, d1, d2, vl, vu, il, iu, abstol, m, wa1, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevx(V,A)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[4 - 1] = ulpinv;
                        result[5 - 1] = ulpinv;
                        result[6 - 1] = ulpinv;
                        goto statement_250;
                    }
                }
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                //
                //              Do tests 4 and 5.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt21(n, 0, d3, d4, wa1, d2, z, ldu, work, &result[4 - 1]);
                //
                ntest = 6;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevx("N", "A", n, d3, d4, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevx(N,A)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[6 - 1] = ulpinv;
                        goto statement_250;
                    }
                }
                //
                //              Do test 6.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(wa2[j - 1]), abs(eveigs[j - 1])});
                    temp2 = max(temp2, abs(wa2[j - 1] - eveigs[j - 1]));
                }
                result[6 - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_250:
                //
                ntest = 7;
                for (i = 1; i <= n; i = i + 1) {
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevr("V", "A", n, d1, d2, vl, vu, il, iu, abstol, m, wa1, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevr(V,A)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[7 - 1] = ulpinv;
                        result[8 - 1] = ulpinv;
                        goto statement_320;
                    }
                }
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                //
                //              Do tests 7 and 8.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt21(n, 0, d3, d4, wa1, d2, z, ldu, work, &result[7 - 1]);
                //
                ntest = 9;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevr("N", "A", n, d3, d4, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevr(N,A)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[9 - 1] = ulpinv;
                        goto statement_320;
                    }
                }
                //
                //              Do test 9.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(wa2[j - 1]), abs(eveigs[j - 1])});
                    temp2 = max(temp2, abs(wa2[j - 1] - eveigs[j - 1]));
                }
                result[9 - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_320:
                //
                ntest = 10;
                for (i = 1; i <= n; i = i + 1) {
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevx("V", "I", n, d1, d2, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevx(V,I)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[10 - 1] = ulpinv;
                        result[11 - 1] = ulpinv;
                        result[12 - 1] = ulpinv;
                        goto statement_380;
                    }
                }
                //
                //              Do tests 10 and 11.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt22(n, m2, 0, d3, d4, wa2, d2, z, ldu, work, max((INTEGER)1, m2), &result[10 - 1]);
                //
                ntest = 12;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevx("N", "I", n, d3, d4, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevx(N,I)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[12 - 1] = ulpinv;
                        goto statement_380;
                    }
                }
                //
                //              Do test 12.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                result[12 - 1] = (temp1 + temp2) / max(unfl, ulp * temp3);
            //
            statement_380:
                //
                ntest = 12;
                if (n > 0) {
                    if (il != 1) {
                        vl = wa1[il - 1] - max({half * (wa1[il - 1] - wa1[(il - 1) - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else {
                        vl = wa1[1 - 1] - max({half * (wa1[n - 1] - wa1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                    if (iu != n) {
                        vu = wa1[iu - 1] + max({half * (wa1[(iu + 1) - 1] - wa1[iu - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else {
                        vu = wa1[n - 1] + max({half * (wa1[n - 1] - wa1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                } else {
                    vl = zero;
                    vu = one;
                }
                //
                for (i = 1; i <= n; i = i + 1) {
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevx("V", "V", n, d1, d2, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevx(V,V)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[13 - 1] = ulpinv;
                        result[14 - 1] = ulpinv;
                        result[15 - 1] = ulpinv;
                        goto statement_440;
                    }
                }
                //
                if (m2 == 0 && n > 0) {
                    result[13 - 1] = ulpinv;
                    result[14 - 1] = ulpinv;
                    result[15 - 1] = ulpinv;
                    goto statement_440;
                }
                //
                //              Do tests 13 and 14.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt22(n, m2, 0, d3, d4, wa2, d2, z, ldu, work, max((INTEGER)1, m2), &result[13 - 1]);
                //
                ntest = 15;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevx("N", "V", n, d3, d4, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevx(N,V)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[15 - 1] = ulpinv;
                        goto statement_440;
                    }
                }
                //
                //              Do test 15.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                result[15 - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
            //
            statement_440:
                //
                ntest = 16;
                for (i = 1; i <= n; i = i + 1) {
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevd("V", n, d1, d2, z, ldu, work, lwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevd(V)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[16 - 1] = ulpinv;
                        result[17 - 1] = ulpinv;
                        result[18 - 1] = ulpinv;
                        goto statement_510;
                    }
                }
                //
                //              Do tests 16 and 17.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt21(n, 0, d3, d4, d1, d2, z, ldu, work, &result[16 - 1]);
                //
                ntest = 18;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevd("N", n, d3, d4, z, ldu, work, lwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevd(N)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[18 - 1] = ulpinv;
                        goto statement_510;
                    }
                }
                //
                //              Do test 18.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(eveigs[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(eveigs[j - 1] - d3[j - 1]));
                }
                result[18 - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_510:
                //
                ntest = 19;
                for (i = 1; i <= n; i = i + 1) {
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevr("V", "I", n, d1, d2, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevr(V,I)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[19 - 1] = ulpinv;
                        result[20 - 1] = ulpinv;
                        result[21 - 1] = ulpinv;
                        goto statement_570;
                    }
                }
                //
                //              DO tests 19 and 20.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt22(n, m2, 0, d3, d4, wa2, d2, z, ldu, work, max((INTEGER)1, m2), &result[19 - 1]);
                //
                ntest = 21;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevr("N", "I", n, d3, d4, vl, vu, il, iu, abstol, m3, wa3, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevr(N,I)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[21 - 1] = ulpinv;
                        goto statement_570;
                    }
                }
                //
                //              Do test 21.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                result[21 - 1] = (temp1 + temp2) / max(unfl, ulp * temp3);
            //
            statement_570:
                //
                ntest = 21;
                if (n > 0) {
                    if (il != 1) {
                        vl = wa1[il - 1] - max({half * (wa1[il - 1] - wa1[(il - 1) - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else {
                        vl = wa1[1 - 1] - max({half * (wa1[n - 1] - wa1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                    if (iu != n) {
                        vu = wa1[iu - 1] + max({half * (wa1[(iu + 1) - 1] - wa1[iu - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else {
                        vu = wa1[n - 1] + max({half * (wa1[n - 1] - wa1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                } else {
                    vl = zero;
                    vu = one;
                }
                //
                for (i = 1; i <= n; i = i + 1) {
                    d1[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d2[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevr("V", "V", n, d1, d2, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevr(V,V)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[22 - 1] = ulpinv;
                        result[23 - 1] = ulpinv;
                        result[24 - 1] = ulpinv;
                        goto statement_630;
                    }
                }
                //
                if (m2 == 0 && n > 0) {
                    result[22 - 1] = ulpinv;
                    result[23 - 1] = ulpinv;
                    result[24 - 1] = ulpinv;
                    goto statement_630;
                }
                //
                //              Do tests 22 and 23.
                //
                for (i = 1; i <= n; i = i + 1) {
                    d3[i - 1] = a[(i - 1) + (i - 1) * lda];
                }
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstt22(n, m2, 0, d3, d4, wa2, d2, z, ldu, work, max((INTEGER)1, m2), &result[22 - 1]);
                //
                ntest = 24;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d4[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                }
                Rstevr("N", "V", n, d3, d4, vl, vu, il, iu, abstol, m3, wa3, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rstevr(N,V)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[24 - 1] = ulpinv;
                        goto statement_630;
                    }
                }
                //
                //              Do test 24.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                result[24 - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
            //
            statement_630:;
                //
            } else {
                //
                for (i = 1; i <= 24; i = i + 1) {
                    result[i - 1] = zero;
                }
                ntest = 24;
            }
            //
            //           Perform remaining tests storing upper or lower triangular
            //           part of matrix.
            //
            for (iuplo = 0; iuplo <= 1; iuplo = iuplo + 1) {
                if (iuplo == 0) {
                    uplo = 'L';
                } else {
                    uplo = 'U';
                }
                //
                //              4)      Call Rsyev and Rsyevx.
                //
                Rlacpy(" ", n, n, a, lda, v, ldu);
                //
                ntest++;
                Rsyev("V", &uplo, n, a, ldu, d1, work, lwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyev(V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyev(V,N)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_660;
                    }
                }
                //
                //              Do tests 25 and 26 (or +54)
                //
                Rsyt21(1, &uplo, n, 0, v, ldu, d1, d2, a, ldu, z, ldu, tau, work, &result[ntest - 1]);
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest += 2;
                Rsyev_2stage("N", &uplo, n, a, ldu, d3, work, lwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyev_2stage(N,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyev_2stage(N,N)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_660;
                    }
                }
                //
                //              Do test 27 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(d1[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(d1[j - 1] - d3[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_660:
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest++;
                //
                if (n > 0) {
                    temp3 = max(abs(d1[1 - 1]), abs(d1[n - 1]));
                    if (il != 1) {
                        vl = d1[il - 1] - max({half * (d1[il - 1] - d1[(il - 1) - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else if (n > 0) {
                        vl = d1[1 - 1] - max({half * (d1[n - 1] - d1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                    if (iu != n) {
                        vu = d1[iu - 1] + max({half * (d1[(iu + 1) - 1] - d1[iu - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else if (n > 0) {
                        vu = d1[n - 1] + max({half * (d1[n - 1] - d1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                } else {
                    temp3 = zero;
                    vl = zero;
                    vu = one;
                }
                //
                Rsyevx("V", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m, wa1, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevx(V,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevx(V,A,N)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_680;
                    }
                }
                //
                //              Do tests 28 and 29 (or +54)
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                Rsyt21(1, &uplo, n, 0, a, ldu, d1, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                Rsyevx_2stage("N", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevx_2stage(N,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevx_2stage(N,A,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_680;
                    }
                }
                //
                //              Do test 30 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(wa1[j - 1]), abs(wa2[j - 1])});
                    temp2 = max(temp2, abs(wa1[j - 1] - wa2[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_680:
                //
                ntest++;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevx("V", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevx(V,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevx(V,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_690;
                    }
                }
                //
                //              Do tests 31 and 32 (or +54)
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevx_2stage("N", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevx_2stage(N,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevx_2stage(N,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_690;
                    }
                }
                //
                //              Do test 33 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                result[ntest - 1] = (temp1 + temp2) / max(unfl, ulp * temp3);
            statement_690:
                //
                ntest++;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevx("V", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevx(V,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevx(V,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_700;
                    }
                }
                //
                //              Do tests 34 and 35 (or +54)
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevx_2stage("N", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevx_2stage(N,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevx_2stage(N,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_700;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_700;
                }
                //
                //              Do test 36 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
            //
            statement_700:
                //
                //              5)      Call Rspev and Rspevx.
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                //              Load array WORK with the upper or lower triangular
                //              part of the matrix in packed form.
                //
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                ntest++;
                Rspev("V", &uplo, n, work, d1, z, ldu, v, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspev(V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspev(V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_800;
                    }
                }
                //
                //              Do tests 37 and 38 (or +54)
                //
                Rsyt21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                ntest += 2;
                Rspev("N", &uplo, n, work, d3, z, ldu, v, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspev(N,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspev(N,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_800;
                    }
                }
                //
                //              Do test 39 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(d1[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(d1[j - 1] - d3[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            //              Load array WORK with the upper or lower triangular part
            //              of the matrix in packed form.
            //
            statement_800:
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                ntest++;
                //
                if (n > 0) {
                    temp3 = max(abs(d1[1 - 1]), abs(d1[n - 1]));
                    if (il != 1) {
                        vl = d1[il - 1] - max({half * (d1[il - 1] - d1[(il - 1) - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else if (n > 0) {
                        vl = d1[1 - 1] - max({half * (d1[n - 1] - d1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                    if (iu != n) {
                        vu = d1[iu - 1] + max({half * (d1[(iu + 1) - 1] - d1[iu - 1]), ten * ulp * temp3, ten * rtunfl});
                    } else if (n > 0) {
                        vu = d1[n - 1] + max({half * (d1[n - 1] - d1[1 - 1]), ten * ulp * temp3, ten * rtunfl});
                    }
                } else {
                    temp3 = zero;
                    vl = zero;
                    vu = one;
                }
                //
                Rspevx("V", "A", &uplo, n, work, vl, vu, il, iu, abstol, m, wa1, z, ldu, v, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevx(V,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevx(V,A,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_900;
                    }
                }
                //
                //              Do tests 40 and 41 (or +54)
                //
                Rsyt21(1, &uplo, n, 0, a, ldu, wa1, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                //
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                Rspevx("N", "A", &uplo, n, work, vl, vu, il, iu, abstol, m2, wa2, z, ldu, v, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevx(N,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevx(N,A,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_900;
                    }
                }
                //
                //              Do test 42 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(wa1[j - 1]), abs(wa2[j - 1])});
                    temp2 = max(temp2, abs(wa1[j - 1] - wa2[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_900:
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                ntest++;
                //
                Rspevx("V", "I", &uplo, n, work, vl, vu, il, iu, abstol, m2, wa2, z, ldu, v, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevx(V,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevx(V,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_990;
                    }
                }
                //
                //              Do tests 43 and 44 (or +54)
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                //
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                Rspevx("N", "I", &uplo, n, work, vl, vu, il, iu, abstol, m3, wa3, z, ldu, v, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevx(N,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevx(N,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_990;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_990;
                }
                //
                //              Do test 45 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
            //
            statement_990:
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                ntest++;
                //
                Rspevx("V", "V", &uplo, n, work, vl, vu, il, iu, abstol, m2, wa2, z, ldu, v, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevx(V,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevx(V,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1080;
                    }
                }
                //
                //              Do tests 46 and 47 (or +54)
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                //
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                Rspevx("N", "V", &uplo, n, work, vl, vu, il, iu, abstol, m3, wa3, z, ldu, v, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevx(N,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevx(N,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1080;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_1080;
                }
                //
                //              Do test 48 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
            //
            statement_1080:
                //
                //              6)      Call Rsbev and Rsbevx.
                //
                if (jtype <= 7) {
                    kd = 1;
                } else if (jtype >= 8 && jtype <= 15) {
                    kd = max(n - 1, 0);
                } else {
                    kd = ihbw;
                }
                //
                //              Load array V with the upper or lower triangular part
                //              of the matrix in band form.
                //
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                ntest++;
                Rsbev("V", &uplo, n, kd, v, ldu, d1, z, ldu, work, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbev(V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbev(V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1180;
                    }
                }
                //
                //              Do tests 49 and 50 (or ... )
                //
                Rsyt21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                ntest += 2;
                Rsbev_2stage("N", &uplo, n, kd, v, ldu, d3, z, ldu, work, lwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbev_2stage(N,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbev_2stage(N,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1180;
                    }
                }
                //
                //              Do test 51 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(d1[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(d1[j - 1] - d3[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            //              Load array V with the upper or lower triangular part
            //              of the matrix in band form.
            //
            statement_1180:
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                ntest++;
                Rsbevx("V", "A", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m, wa2, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevx(V,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevx(V,A,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1280;
                    }
                }
                //
                //              Do tests 52 and 53 (or +54)
                //
                Rsyt21(1, &uplo, n, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                //
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                Rsbevx_2stage("N", "A", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevx_2stage(N,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevx_2stage(N,A,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1280;
                    }
                }
                //
                //              Do test 54 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(wa2[j - 1]), abs(wa3[j - 1])});
                    temp2 = max(temp2, abs(wa2[j - 1] - wa3[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_1280:
                ntest++;
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                Rsbevx("V", "I", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevx(V,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevx(V,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1370;
                    }
                }
                //
                //              Do tests 55 and 56 (or +54)
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                //
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                Rsbevx_2stage("N", "I", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevx_2stage(N,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevx_2stage(N,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1370;
                    }
                }
                //
                //              Do test 57 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
            //
            statement_1370:
                ntest++;
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                Rsbevx("V", "V", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevx(V,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevx(V,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1460;
                    }
                }
                //
                //              Do tests 58 and 59 (or +54)
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                //
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                Rsbevx_2stage("N", "V", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, lwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevx_2stage(N,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevx_2stage(N,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1460;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_1460;
                }
                //
                //              Do test 60 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
            //
            statement_1460:
                //
                //              7)      Call Rsyevd
                //
                Rlacpy(" ", n, n, a, lda, v, ldu);
                //
                ntest++;
                Rsyevd("V", &uplo, n, a, ldu, d1, work, lwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevd(V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevd(V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1480;
                    }
                }
                //
                //              Do tests 61 and 62 (or +54)
                //
                Rsyt21(1, &uplo, n, 0, v, ldu, d1, d2, a, ldu, z, ldu, tau, work, &result[ntest - 1]);
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest += 2;
                Rsyevd_2stage("N", &uplo, n, a, ldu, d3, work, lwork, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevd_2stage(N,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevd_2stage(N,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1480;
                    }
                }
                //
                //              Do test 63 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(d1[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(d1[j - 1] - d3[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_1480:
                //
                //              8)      Call Rspevd.
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                //              Load array WORK with the upper or lower triangular
                //              part of the matrix in packed form.
                //
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                ntest++;
                Rspevd("V", &uplo, n, work, d1, z, ldu, &work[indx - 1], lwedc - indx + 1, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevd(V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevd(V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1580;
                    }
                }
                //
                //              Do tests 64 and 65 (or +54)
                //
                Rsyt21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                if (iuplo == 1) {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            //
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                } else {
                    indx = 1;
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            work[indx - 1] = a[(i - 1) + (j - 1) * lda];
                            indx++;
                        }
                    }
                }
                //
                ntest += 2;
                Rspevd("N", &uplo, n, work, d3, z, ldu, &work[indx - 1], lwedc - indx + 1, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rspevd(N,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rspevd(N,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1580;
                    }
                }
                //
                //              Do test 66 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(d1[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(d1[j - 1] - d3[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            statement_1580:
                //
                //              9)      Call Rsbevd.
                //
                if (jtype <= 7) {
                    kd = 1;
                } else if (jtype >= 8 && jtype <= 15) {
                    kd = max(n - 1, 0);
                } else {
                    kd = ihbw;
                }
                //
                //              Load array V with the upper or lower triangular part
                //              of the matrix in band form.
                //
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                ntest++;
                Rsbevd("V", &uplo, n, kd, v, ldu, d1, z, ldu, work, lwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevd(V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevd(V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1680;
                    }
                }
                //
                //              Do tests 67 and 68 (or +54)
                //
                Rsyt21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                if (iuplo == 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = max((INTEGER)1, j - kd); i <= j; i = i + 1) {
                            v[((kd + 1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= min(n, j + kd); i = i + 1) {
                            v[((1 + i - j) - 1) + (j - 1) * ldv] = a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
                //
                ntest += 2;
                Rsbevd_2stage("N", &uplo, n, kd, v, ldu, d3, z, ldu, work, lwork, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsbevd_2stage(N,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsbevd_2stage(N,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1680;
                    }
                }
                //
                //              Do test 69 (or +54)
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(d1[j - 1]), abs(d3[j - 1])});
                    temp2 = max(temp2, abs(d1[j - 1] - d3[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_1680:
                //
                Rlacpy(" ", n, n, a, lda, v, ldu);
                ntest++;
                Rsyevr("V", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m, wa1, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevr(V,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevr(V,A,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1700;
                    }
                }
                //
                //              Do tests 70 and 71 (or ... )
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                Rsyt21(1, &uplo, n, 0, a, ldu, wa1, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                Rsyevr_2stage("N", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevr_2stage(N,A,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevr_2stage(N,A,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1700;
                    }
                }
                //
                //              Do test 72 (or ... )
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, abs(wa1[j - 1]), abs(wa2[j - 1])});
                    temp2 = max(temp2, abs(wa1[j - 1] - wa2[j - 1]));
                }
                result[ntest - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            statement_1700:
                //
                ntest++;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevr("V", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevr(V,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevr(V,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1710;
                    }
                }
                //
                //              Do tests 73 and 74 (or +54)
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevr_2stage("N", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevr_2stage(N,I,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevr_2stage(N,I,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1710;
                    }
                }
                //
                //              Do test 75 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                result[ntest - 1] = (temp1 + temp2) / max(unfl, ulp * temp3);
            statement_1710:
                //
                ntest++;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevr("V", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevr(V,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevr(V,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_700;
                    }
                }
                //
                //              Do tests 76 and 77 (or +54)
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
                Rsyt22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, &result[ntest - 1]);
                //
                ntest += 2;
                Rlacpy(" ", n, n, v, ldu, a, lda);
                Rsyevr_2stage("N", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, iwork, work, lwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Rsyevr_2stage(N,V,U)", iinfo, n, jtype, ioldsd;
                    else
                        write(nounit, format_9999), "Rsyevr_2stage(N,V,L)", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_700;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_700;
                }
                //
                //              Do test 78 (or +54)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, temp3 * ulp);
                //
                Rlacpy(" ", n, n, v, ldu, a, lda);
                //
            }
            //
            //           End of Loop -- Check for RESULT(j) > THRESH
            //
            ntestt += ntest;
            //
            Rlafts("DST", n, n, jtype, ntest, result, ioldsd, thresh, nounit, nerrs);
        //
        statement_1730:;
        }
    }
    //
    //     Summary
    //
    Alasvm("DST", nounit, nerrs, ntestt, 0);
    //
    //     End of Rdrvst2stg
    //
}

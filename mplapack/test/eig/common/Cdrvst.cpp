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

void Cdrvst(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, COMPLEX *a, INTEGER const lda, REAL *d1, REAL *d2, REAL *d3, REAL *wa1, REAL *wa2, REAL *wa3, COMPLEX *u, INTEGER const ldu, COMPLEX *v, COMPLEX *tau, COMPLEX *z, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER *iwork, INTEGER const liwork, REAL *result, INTEGER &info) {

    INTEGER ldv = ldu;
    INTEGER ldz = ldu;

    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 18;
    INTEGER ktype[18] = {1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9};
    INTEGER kmagn[18] = {1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3};
    INTEGER kmode[18] = {0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4};
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
    const REAL two = 2.0e+0;
    INTEGER lgn = 0;
    INTEGER lwedc = 0;
    INTEGER lrwedc = 0;
    INTEGER liwedc = 0;
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
    INTEGER jcol = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER idumma[1];
    const REAL zero = 0.0;
    INTEGER ihbw = 0;
    INTEGER idiag = 0;
    INTEGER irow = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    REAL abstol = 0.0;
    INTEGER il = 0;
    INTEGER iu = 0;
    INTEGER itemp = 0;
    INTEGER iuplo = 0;
    char uplo;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    REAL temp3 = 0.0;
    const REAL half = one / two;
    const REAL ten = 10.0;
    REAL vl = 0.0;
    REAL vu = 0.0;
    INTEGER m = 0;
    INTEGER m2 = 0;
    INTEGER m3 = 0;
    INTEGER indx = 0;
    INTEGER indwrk = 0;
    INTEGER kd = 0;
    static const char *format_9998 = "(' Cdrvst: ',a,' returned INFO=',i6,/,9x,'N=',i6,', KD=',i6,', JTYPE=',"
                                     "i6,', ISEED=(',3(i5,','),i5,')')";
    static const char *format_9999 = "(' Cdrvst: ',a,' returned INFO=',i6,/,9x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',3(i5,','),i5,')')";
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
    } else if (2 * max((INTEGER)2, nmax) * max((INTEGER)2, nmax) > lwork) {
        info = -22;
    }
    //
    if (info != 0) {
        Mxerbla("Cdrvst", -info);
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
            if ((INTEGER)pow((double)2, (double)lgn) < n) {
                lgn++;
            }
            if ((INTEGER)pow((double)2, (double)lgn) < n) {
                lgn++;
            }
            lwedc = max((INTEGER)2 * n + n * n, 2 * n * n);
            lrwedc = 1 + 4 * n + 2 * n * lgn + 3 * n * n;
            liwedc = 3 + 5 * n;
        } else {
            lwedc = 2;
            lrwedc = 8;
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
            if (!dotype[jtype - 1]) {
                goto statement_1210;
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
            //           =5         random log   Hermitian, w/ eigenvalues
            //           =6         random       (none)
            //           =7                      random diagonal
            //           =8                      random Hermitian
            //           =9                      band Hermitian, w/ eigenvalues
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
            Claset("Full", lda, n, czero, czero, a, lda);
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
                //              Hermitian banded, eigenvalues specified
                //
                ihbw = castINTEGER((n - 1) * Rlarnd(1, iseed3));
                Clatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, ihbw, ihbw, "Z", u, ldu, work, iinfo);
                //
                //              Store as dense matrix for most routines.
                //
                Claset("Full", lda, n, czero, czero, a, lda);
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
                write(nounit, format_9999), "Generator", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
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
            //           Perform tests storing upper or lower triangular
            //           part of matrix.
            //
            for (iuplo = 0; iuplo <= 1; iuplo = iuplo + 1) {
                if (iuplo == 0) {
                    uplo = 'L';
                } else {
                    uplo = 'U';
                }
                //
                //              Call Cheevd and CHEEVX.
                //
                Clacpy(" ", n, n, a, lda, v, ldu);
                //
                ntest++;
                Cheevd("V", &uplo, n, a, ldu, d1, work, lwedc, rwork, lrwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevd(V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevd(V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_130;
                    }
                }
                //
                //              Do tests 1 and 2.
                //
                Chet21(1, &uplo, n, 0, v, ldu, d1, d2, a, ldu, z, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest += 2;
                Cheevd("N", &uplo, n, a, ldu, d3, work, lwedc, rwork, lrwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevd(N,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevd(N,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_130;
                    }
                }
                //
                //              Do test 3.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d1[j - 1] - d3[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            statement_130:
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest++;
                //
                if (n > 0) {
                    temp3 = max(abs(d1[1 - 1]), abs(d1[n - 1]));
                    if (il != 1) {
                        vl = d1[il - 1] - max({REAL(half * (d1[il - 1] - d1[(il - 1) - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    } else if (n > 0) {
                        vl = d1[1 - 1] - max({REAL(half * (d1[n - 1] - d1[1 - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    }
                    if (iu != n) {
                        vu = d1[iu - 1] + max({REAL(half * (d1[(iu + 1) - 1] - d1[iu - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    } else if (n > 0) {
                        vu = d1[n - 1] + max({REAL(half * (d1[n - 1] - d1[1 - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    }
                } else {
                    temp3 = zero;
                    vl = zero;
                    vu = one;
                }
                //
                Cheevx("V", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m, wa1, z, ldu, work, lwork, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevx(V,A,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevx(V,A,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_150;
                    }
                }
                //
                //              Do tests 4 and 5.
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                Chet21(1, &uplo, n, 0, a, ldu, wa1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                ntest += 2;
                Cheevx("N", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, lwork, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevx(N,A,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevx(N,A,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_150;
                    }
                }
                //
                //              Do test 6.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(wa1[j - 1])), REAL(abs(wa2[j - 1]))});
                    temp2 = max(temp2, REAL(abs(wa1[j - 1] - wa2[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            statement_150:
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest++;
                //
                Cheevx("V", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, lwork, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevx(V,I,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevx(V,I,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_160;
                    }
                }
                //
                //              Do tests 7 and 8.
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                ntest += 2;
                //
                Cheevx("N", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, lwork, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevx(N,I,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevx(N,I,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_160;
                    }
                }
                //
                //              Do test 9.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
            //
            statement_160:
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest++;
                //
                Cheevx("V", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, lwork, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevx(V,V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevx(V,V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_170;
                    }
                }
                //
                //              Do tests 10 and 11.
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                ntest += 2;
                //
                Cheevx("N", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, lwork, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevx(N,V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevx(N,V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_170;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_170;
                }
                //
                //              Do test 12.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
            //
            statement_170:
                //
                //              Call Chpevd and CHPEVX.
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
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
                indwrk = n * (n + 1) / 2 + 1;
                Chpevd("V", &uplo, n, work, d1, z, ldu, &work[indwrk - 1], lwedc, rwork, lrwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevd(V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevd(V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_270;
                    }
                }
                //
                //              Do tests 13 and 14.
                //
                Chet21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                indwrk = n * (n + 1) / 2 + 1;
                Chpevd("N", &uplo, n, work, d3, z, ldu, &work[indwrk - 1], lwedc, rwork, lrwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevd(N,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevd(N,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_270;
                    }
                }
                //
                //              Do test 15.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d1[j - 1] - d3[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            //              Load array WORK with the upper or lower triangular part
            //              of the matrix in packed form.
            //
            statement_270:
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
                        vl = d1[il - 1] - max({REAL(half * (d1[il - 1] - d1[(il - 1) - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    } else if (n > 0) {
                        vl = d1[1 - 1] - max({REAL(half * (d1[n - 1] - d1[1 - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    }
                    if (iu != n) {
                        vu = d1[iu - 1] + max({REAL(half * (d1[(iu + 1) - 1] - d1[iu - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    } else if (n > 0) {
                        vu = d1[n - 1] + max({REAL(half * (d1[n - 1] - d1[1 - 1])), REAL(ten * ulp * temp3), REAL(ten * rtunfl)});
                    }
                } else {
                    temp3 = zero;
                    vl = zero;
                    vu = one;
                }
                //
                Chpevx("V", "A", &uplo, n, work, vl, vu, il, iu, abstol, m, wa1, z, ldu, v, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevx(V,A,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevx(V,A,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_370;
                    }
                }
                //
                //              Do tests 16 and 17.
                //
                Chet21(1, &uplo, n, 0, a, ldu, wa1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                Chpevx("N", "A", &uplo, n, work, vl, vu, il, iu, abstol, m2, wa2, z, ldu, v, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevx(N,A,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevx(N,A,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_370;
                    }
                }
                //
                //              Do test 18.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(wa1[j - 1])), REAL(abs(wa2[j - 1]))});
                    temp2 = max(temp2, REAL(abs(wa1[j - 1] - wa2[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            statement_370:
                ntest++;
                if (i & uplo == 1) {
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
                Chpevx("V", "I", &uplo, n, work, vl, vu, il, iu, abstol, m2, wa2, z, ldu, v, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevx(V,I,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevx(V,I,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_460;
                    }
                }
                //
                //              Do tests 19 and 20.
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                ntest += 2;
                //
                if (i & uplo == 1) {
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
                Chpevx("N", "I", &uplo, n, work, vl, vu, il, iu, abstol, m3, wa3, z, ldu, v, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevx(N,I,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevx(N,I,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_460;
                    }
                }
                //
                //              Do test 21.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
            //
            statement_460:
                ntest++;
                if (i & uplo == 1) {
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
                Chpevx("V", "V", &uplo, n, work, vl, vu, il, iu, abstol, m2, wa2, z, ldu, v, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevx(V,V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevx(V,V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_550;
                    }
                }
                //
                //              Do tests 22 and 23.
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                Chpevx("N", "V", &uplo, n, work, vl, vu, il, iu, abstol, m3, wa3, z, ldu, v, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpevx(N,V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpevx(N,V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_550;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_550;
                }
                //
                //              Do test 24.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
            //
            statement_550:
                //
                //              Call Chbevd and CHBEVX.
                //
                if (jtype <= 7) {
                    kd = 0;
                } else if (jtype >= 8 && jtype <= 15) {
                    kd = max(n - 1, (INTEGER)0);
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
                Chbevd("V", &uplo, n, kd, v, ldu, d1, z, ldu, work, lwedc, rwork, lrwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbevd(V,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbevd(V,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_650;
                    }
                }
                //
                //              Do tests 25 and 26.
                //
                Chet21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                Chbevd("N", &uplo, n, kd, v, ldu, d3, z, ldu, work, lwedc, rwork, lrwedc, iwork, liwedc, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbevd(N,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbevd(N,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_650;
                    }
                }
                //
                //              Do test 27.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d1[j - 1] - d3[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            //              Load array V with the upper or lower triangular part
            //              of the matrix in band form.
            //
            statement_650:
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
                Chbevx("V", "A", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m, wa1, z, ldu, work, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chbevx(V,A,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chbevx(V,A,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_750;
                    }
                }
                //
                //              Do tests 28 and 29.
                //
                Chet21(1, &uplo, n, 0, a, ldu, wa1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                Chbevx("N", "A", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbevx(N,A,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbevx(N,A,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_750;
                    }
                }
                //
                //              Do test 30.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(wa1[j - 1])), REAL(abs(wa2[j - 1]))});
                    temp2 = max(temp2, REAL(abs(wa1[j - 1] - wa2[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            //              Load array V with the upper or lower triangular part
            //              of the matrix in band form.
            //
            statement_750:
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
                Chbevx("V", "I", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbevx(V,I,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbevx(V,I,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_840;
                    }
                }
                //
                //              Do tests 31 and 32.
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                Chbevx("N", "I", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbevx(N,I,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbevx(N,I,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_840;
                    }
                }
                //
                //              Do test 33.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
            //
            //              Load array V with the upper or lower triangular part
            //              of the matrix in band form.
            //
            statement_840:
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
                Chbevx("V", "V", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, work, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbevx(V,V,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbevx(V,V,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_930;
                    }
                }
                //
                //              Do tests 34 and 35.
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                Chbevx("N", "V", &uplo, n, kd, v, ldu, u, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, work, rwork, iwork, &iwork[(5 * n + 1) - 1], iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbevx(N,V,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbevx(N,V,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_930;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_930;
                }
                //
                //              Do test 36.
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
            //
            statement_930:
                //
                //              Call Cheev
                //
                Clacpy(" ", n, n, a, lda, v, ldu);
                //
                ntest++;
                Cheev("V", &uplo, n, a, ldu, d1, work, lwork, rwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheev(V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheev(V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_950;
                    }
                }
                //
                //              Do tests 37 and 38
                //
                Chet21(1, &uplo, n, 0, v, ldu, d1, d2, a, ldu, z, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                ntest += 2;
                Cheev("N", &uplo, n, a, ldu, d3, work, lwork, rwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheev(N,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheev(N,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_950;
                    }
                }
                //
                //              Do test 39
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d1[j - 1] - d3[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            statement_950:
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                //              Call Chpev
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
                indwrk = n * (n + 1) / 2 + 1;
                Chpev("V", &uplo, n, work, d1, z, ldu, &work[indwrk - 1], rwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpev(V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpev(V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1050;
                    }
                }
                //
                //              Do tests 40 and 41.
                //
                Chet21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                indwrk = n * (n + 1) / 2 + 1;
                Chpev("N", &uplo, n, work, d3, z, ldu, &work[indwrk - 1], rwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Chpev(N,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Chpev(N,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1050;
                    }
                }
                //
                //              Do test 42
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d1[j - 1] - d3[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            statement_1050:
                //
                //              Call Chbev
                //
                if (jtype <= 7) {
                    kd = 0;
                } else if (jtype >= 8 && jtype <= 15) {
                    kd = max(n - 1, (INTEGER)0);
                } else {
                    kd = ihbw;
                }
                //
                //              Load array V with the upper or lower triangular part
                //              of the matrix in band form.
                //
                if (i & uplo == 1) {
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
                Chbev("V", &uplo, n, kd, v, ldu, d1, z, ldu, work, rwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbev(V,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbev(V,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1140;
                    }
                }
                //
                //              Do tests 43 and 44.
                //
                Chet21(1, &uplo, n, 0, a, lda, d1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
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
                Chbev("N", &uplo, n, kd, v, ldu, d3, z, ldu, work, rwork, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9998), "Chbev(N,U)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9998), "Chbev(N,L)", iinfo, n, kd, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1140;
                    }
                }
            //
            statement_1140:
                //
                //              Do test 45.
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(d1[j - 1])), REAL(abs(d3[j - 1]))});
                    temp2 = max(temp2, REAL(abs(d1[j - 1] - d3[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
                //
                Clacpy(" ", n, n, a, lda, v, ldu);
                ntest++;
                Cheevr("V", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m, wa1, z, ldu, iwork, work, lwork, rwork, lrwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevr(V,A,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevr(V,A,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1170;
                    }
                }
                //
                //              Do tests 45 and 46 (or ... )
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                Chet21(1, &uplo, n, 0, a, ldu, wa1, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                ntest += 2;
                Cheevr("N", "A", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, rwork, lrwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevr(N,A,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevr(N,A,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1170;
                    }
                }
                //
                //              Do test 47 (or ... )
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max({temp1, REAL(abs(wa1[j - 1])), REAL(abs(wa2[j - 1]))});
                    temp2 = max(temp2, REAL(abs(wa1[j - 1] - wa2[j - 1])));
                }
                result[ntest - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            statement_1170:
                //
                ntest++;
                Clacpy(" ", n, n, v, ldu, a, lda);
                Cheevr("V", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, rwork, lrwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevr(V,I,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevr(V,I,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
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
                //              Do tests 48 and 49 (or +??)
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                ntest += 2;
                Clacpy(" ", n, n, v, ldu, a, lda);
                Cheevr("N", "I", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, iwork, work, lwork, rwork, lrwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevr(N,I,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevr(N,I,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1180;
                    }
                }
                //
                //              Do test 50 (or +??)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(ulp * temp3));
            statement_1180:
                //
                ntest++;
                Clacpy(" ", n, n, v, ldu, a, lda);
                Cheevr("V", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m2, wa2, z, ldu, iwork, work, lwork, rwork, lrwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevr(V,V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevr(V,V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        result[(ntest + 1) - 1] = ulpinv;
                        result[(ntest + 2) - 1] = ulpinv;
                        goto statement_1190;
                    }
                }
                //
                //              Do tests 51 and 52 (or +??)
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
                //
                Chet22(1, &uplo, n, m2, 0, a, ldu, wa2, d2, z, ldu, v, ldu, tau, work, rwork, &result[ntest - 1]);
                //
                ntest += 2;
                Clacpy(" ", n, n, v, ldu, a, lda);
                Cheevr("N", "V", &uplo, n, a, ldu, vl, vu, il, iu, abstol, m3, wa3, z, ldu, iwork, work, lwork, rwork, lrwork, &iwork[(2 * n + 1) - 1], liwork - 2 * n, iinfo);
                if (iinfo != 0) {
                    if (Mlsame(&uplo, "U"))
                        write(nounit, format_9999), "Cheevr(N,V,U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    else
                        write(nounit, format_9999), "Cheevr(N,V,L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[ntest - 1] = ulpinv;
                        goto statement_1190;
                    }
                }
                //
                if (m3 == 0 && n > 0) {
                    result[ntest - 1] = ulpinv;
                    goto statement_1190;
                }
                //
                //              Do test 52 (or +??)
                //
                temp1 = Rsxt1(1, wa2, m2, wa3, m3, abstol, ulp, unfl);
                temp2 = Rsxt1(1, wa3, m3, wa2, m2, abstol, ulp, unfl);
                if (n > 0) {
                    temp3 = max(abs(wa1[1 - 1]), abs(wa1[n - 1]));
                } else {
                    temp3 = zero;
                }
                result[ntest - 1] = (temp1 + temp2) / max(unfl, REAL(temp3 * ulp));
                //
                Clacpy(" ", n, n, v, ldu, a, lda);
            //
            //              Load array V with the upper or lower triangular part
            //              of the matrix in band form.
            //
            statement_1190:;
                //
            }
            //
            //           End of Loop -- Check for RESULT(j) > THRESH
            //
            ntestt += ntest;
            Rlafts("ZST", n, n, jtype, ntest, result, ioldsd, thresh, nounit, nerrs);
        //
        statement_1210:;
        }
    }
    //
    //     Summary
    //
    Alasvm("ZST", nounit, nerrs, ntestt, 0);
    //
    //     End of Cdrvst
    //
}

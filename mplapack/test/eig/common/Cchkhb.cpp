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

void Cchkhb(INTEGER const nsizes, INTEGER *nn, INTEGER const nwdths, INTEGER *kk, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, COMPLEX *a, INTEGER const lda, REAL *sd, REAL *se, COMPLEX *u, INTEGER const ldu, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result, INTEGER &info) {
    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 15;
    INTEGER ktype[15] = {1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8};
    INTEGER kmagn[15] = {1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3};
    INTEGER kmode[15] = {0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0};
    INTEGER ntestt = 0;
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    bool badnnb = false;
    INTEGER kmax = 0;
    REAL unfl = 0.0;
    const REAL one = 1.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    REAL ulpinv = 0.0;
    REAL rtunfl = 0.0;
    REAL rtovfl = 0.0;
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    REAL aninv = 0.0;
    INTEGER jwidth = 0;
    INTEGER k = 0;
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
    INTEGER jcol = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER idumma[1];
    const REAL zero = 0.0;
    INTEGER i = 0;
    REAL temp1 = 0.0;
    const REAL two = 2.0e+0;
    const REAL half = one / two;
    INTEGER jc = 0;
    INTEGER jr = 0;
    char buf[1024];
    static const char *format_9999 = "(' Cchkhb: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
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
    //     Check for errors
    //
    ntestt = 0;
    info = 0;
    //
    //     Important constants
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
    badnnb = false;
    kmax = 0;
    for (j = 1; j <= nsizes; j = j + 1) {
        kmax = max(kmax, kk[j - 1]);
        if (kk[j - 1] < 0) {
            badnnb = true;
        }
    }
    kmax = min(nmax - 1, kmax);
    //
    //     Check for errors
    //
    if (nsizes < 0) {
        info = -1;
    } else if (badnn) {
        info = -2;
    } else if (nwdths < 0) {
        info = -3;
    } else if (badnnb) {
        info = -4;
    } else if (ntypes < 0) {
        info = -5;
    } else if (lda < kmax + 1) {
        info = -11;
    } else if (ldu < nmax) {
        info = -15;
    } else if ((max(lda, nmax) + 1) * nmax > lwork) {
        info = -17;
    }
    //
    if (info != 0) {
        Mxerbla("Cchkhb", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (nsizes == 0 || ntypes == 0 || nwdths == 0) {
        return;
    }
    //
    //     More Important constants
    //
    unfl = Rlamch("Safe minimum");
    ovfl = one / unfl;
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    ulpinv = one / ulp;
    rtunfl = sqrt(unfl);
    rtovfl = sqrt(ovfl);
    //
    //     Loop over sizes, types
    //
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        aninv = one / castREAL(max((INTEGER)1, n));
        //
        for (jwidth = 1; jwidth <= nwdths; jwidth = jwidth + 1) {
            k = kk[jwidth - 1];
            if (k > n) {
                goto statement_180;
            }
            k = max({(INTEGER)0, min(n - 1, k)});
            //
            if (nsizes != 1) {
                mtypes = min(maxtyp, ntypes);
            } else {
                mtypes = min(maxtyp + 1, ntypes);
            }
            //
            for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
                if (!dotype[jtype - 1]) {
                    goto statement_170;
                }
                nmats++;
                ntest = 0;
                //
                for (j = 1; j <= 4; j = j + 1) {
                    ioldsd[j - 1] = iseed[j - 1];
                }
                //
                //              Compute "A".
                //              Store as "Upper"; later, we will copy to other format.
                //
                //              Control parameters:
                //
                //                  KMAGN  KMODE        KTYPE
                //              =1  O(1)   clustered 1  zero
                //              =2  large  clustered 2  identity
                //              =3  small  exponential  (none)
                //              =4         arithmetic   diagonal, (w/ eigenvalues)
                //              =5         random log   hermitian, w/ eigenvalues
                //              =6         random       (none)
                //              =7                      random diagonal
                //              =8                      random hermitian
                //              =9                      positive definite
                //              =10                     diagonally dominant tridiagonal
                //
                if (mtypes > maxtyp) {
                    goto statement_100;
                }
                //
                itype = ktype[jtype - 1];
                imode = kmode[jtype - 1];
                //
                //              Compute norm
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
                //              Special Matrices -- Identity & Jordan block
                //
                //                 Zero
                //
                if (itype == 1) {
                    iinfo = 0;
                    //
                } else if (itype == 2) {
                    //
                    //                 Identity
                    //
                    for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                        a[((k + 1) - 1) + (jcol - 1) * lda] = anorm;
                    }
                    //
                } else if (itype == 4) {
                    //
                    //                 Diagonal Matrix, [Eigen]values Specified
                    //
                    Clatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, 0, 0, "Q", &a[((k + 1) - 1)], lda, work, iinfo);
                    //
                } else if (itype == 5) {
                    //
                    //                 Hermitian, eigenvalues specified
                    //
                    Clatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, k, k, "Q", a, lda, work, iinfo);
                    //
                } else if (itype == 7) {
                    //
                    //                 Diagonal, random eigenvalues
                    //
                    Clatmr(n, n, "S", iseed, "H", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "Q", &a[((k + 1) - 1)], lda, idumma, iinfo);
                    //
                } else if (itype == 8) {
                    //
                    //                 Hermitian, random eigenvalues
                    //
                    Clatmr(n, n, "S", iseed, "H", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, k, k, zero, anorm, "Q", a, lda, idumma, iinfo);
                    //
                } else if (itype == 9) {
                    //
                    //                 Positive definite, eigenvalues specified.
                    //
                    Clatms(n, n, "S", iseed, "P", rwork, imode, cond, anorm, k, k, "Q", a, lda, &work[(n + 1) - 1], iinfo);
                    //
                } else if (itype == 10) {
                    //
                    //                 Positive definite tridiagonal, eigenvalues specified.
                    //
                    if (n > 1) {
                        k = max((INTEGER)1, k);
                    }
                    Clatms(n, n, "S", iseed, "P", rwork, imode, cond, anorm, 1, 1, "Q", &a[(k - 1)], lda, work, iinfo);
                    for (i = 2; i <= n; i = i + 1) {
                        temp1 = abs(a[(k - 1) + (i - 1) * lda]) / sqrt(abs(a[((k + 1) - 1) + ((i - 1) - 1) * lda] * a[((k + 1) - 1) + (i - 1) * lda]));
                        if (temp1 > half) {
                            a[(k - 1) + (i - 1) * lda] = half * sqrt(abs(a[((k + 1) - 1) + ((i - 1) - 1) * lda] * a[((k + 1) - 1) + (i - 1) * lda]));
                        }
                    }
                    //
                } else {
                    //
                    iinfo = 1;
                }
                //
                if (iinfo != 0) {
                    write(nounit, format_9999), "Generator", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
            //
            statement_100:
                //
                //              Call Chbtrd to compute S and U from upper triangle.
                //
                Clacpy(" ", k + 1, n, a, lda, work, lda);
                //
                ntest = 1;
                Chbtrd("V", "U", n, k, work, lda, sd, se, u, ldu, &work[(lda * n + 1) - 1], iinfo);
                //
                if (iinfo != 0) {
                    write(nounit, format_9999), "Chbtrd(U)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[1 - 1] = ulpinv;
                        goto statement_150;
                    }
                }
                //
                //              Do tests 1 and 2
                //
                Chbt21("Upper", n, k, 1, a, lda, sd, se, u, ldu, work, rwork, &result[1 - 1]);
                //
                //              Convert A from Upper-Triangle-Only storage to
                //              Lower-Triangle-Only storage.
                //
                for (jc = 1; jc <= n; jc = jc + 1) {
                    for (jr = 0; jr <= min(k, n - jc); jr = jr + 1) {
                        a[((jr + 1) - 1) + (jc - 1) * lda] = conj(a[((k + 1 - jr) - 1) + ((jc + jr) - 1) * lda]);
                    }
                }
                for (jc = n + 1 - k; jc <= n; jc = jc + 1) {
                    for (jr = min(k, n - jc) + 1; jr <= k; jr = jr + 1) {
                        a[((jr + 1) - 1) + (jc - 1) * lda] = zero;
                    }
                }
                //
                //              Call Chbtrd to compute S and U from lower triangle
                //
                Clacpy(" ", k + 1, n, a, lda, work, lda);
                //
                ntest = 3;
                Chbtrd("V", "L", n, k, work, lda, sd, se, u, ldu, &work[(lda * n + 1) - 1], iinfo);
                //
                if (iinfo != 0) {
                    write(nounit, format_9999), "Chbtrd(L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[3 - 1] = ulpinv;
                        goto statement_150;
                    }
                }
                ntest = 4;
                //
                //              Do tests 3 and 4
                //
                Chbt21("Lower", n, k, 1, a, lda, sd, se, u, ldu, work, rwork, &result[3 - 1]);
            //
            //              End of Loop -- Check for RESULT(j) > THRESH
            //
            statement_150:
                ntestt += ntest;
                //
                //              Print out tests which fail.
                //
                for (jr = 1; jr <= ntest; jr = jr + 1) {
                    if (result[jr - 1] >= thresh) {
                        //
                        //                    If this is the first test to fail,
                        //                    print a header to the data file.
                        //
                        if (nerrs == 0) {
                            write(nounit, "(/,1x,a3,"
                                          "' -- Complex Hermitian Banded Tridiagonal Reduction Routines'"
                                          ")"),
                                "ZHB";
                            write(nounit, "(' Matrix types (see DCHK23 for details): ')");
                            write(nounit, "(/,' Special Matrices:',/,"
                                          "'  1=Zero matrix.                        ',"
                                          "'  5=Diagonal: clustered entries.',/,"
                                          "'  2=Identity matrix.                    ',"
                                          "'  6=Diagonal: large, evenly spaced.',/,"
                                          "'  3=Diagonal: evenly spaced entries.    ',"
                                          "'  7=Diagonal: small, evenly spaced.',/,"
                                          "'  4=Diagonal: geometr. spaced entries.')");
                            write(nounit, "(' Dense ',a,' Banded Matrices:',/,"
                                          "'  8=Evenly spaced eigenvals.            ',"
                                          "' 12=Small, evenly spaced eigenvals.',/,"
                                          "'  9=Geometrically spaced eigenvals.     ',"
                                          "' 13=Matrix with random O(1) entries.',/,"
                                          "' 10=Clustered eigenvalues.              ',"
                                          "' 14=Matrix with large random entries.',/,"
                                          "' 11=Large, evenly spaced eigenvals.     ',"
                                          "' 15=Matrix with small random entries.')"),
                                "Hermitian";
                            {
                                write_loop wloop(cmn, nounit,
                                                 "(/,' Tests performed:   (S is Tridiag,  U is ',a,',',/,20x,"
                                                 "a,' means ',a,'.',/,' UPLO=''U'':',/,'  1= | A - U S U',a1,"
                                                 "' | / ( |A| n ulp )     ','  2= | I - U U',a1,"
                                                 "' | / ( n ulp )',/,' UPLO=''L'':',/,'  3= | A - U S U',a1,"
                                                 "' | / ( |A| n ulp )     ','  4= | I - U U',a1,"
                                                 "' | / ( n ulp )')");
                                wloop, "unitary", "*", "conjugate transpose";
                                for (j = 1; j <= 4; j = j + 1) {
                                    wloop, "*";
                                }
                            }
                        }
                        nerrs++;
                        sprintnum_short(buf, result[jr - 1]);
                        write(nounit, "(' N=',i5,', K=',i4,', seed=',4(i4,','),' type ',i2,', test(',"
                                      "i2,')=',a)"),
                            n, k, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3], jtype, jr, buf;
                    }
                }
            //
            statement_170:;
            }
        statement_180:;
        }
    }
    //
    //     Summary
    //
    Rlasum("ZHB", nounit, nerrs, ntestt);
    //
    //     End of Cchkhb
    //
}

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

void Rdrves(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *h, REAL *ht, REAL *wr, REAL *wi, REAL *wrt, REAL *wit, REAL *vs, INTEGER const ldvs, REAL *result, REAL *work, INTEGER const nwork, INTEGER *iwork, bool *bwork, INTEGER &info) {

    INTEGER ldh = lda;
    INTEGER ldht = lda;
    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 21;
    INTEGER ktype[21] = {1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9};
    INTEGER kmagn[21] = {1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3};
    INTEGER kmode[21] = {0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1};
    INTEGER kconds[21] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0};
    char path[3];
    char buf[1024];
    INTEGER ntestt = 0;
    INTEGER ntestf = 0;
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL unfl = 0.0;
    const REAL one = 1.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    REAL ulpinv = 0.0;
    REAL rtulp = 0.0;
    REAL rtulpi = 0.0;
    INTEGER nerrs = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ioldsd[4];
    INTEGER itype = 0;
    INTEGER imode = 0;
    REAL anorm = 0.0;
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    INTEGER jcol = 0;
    REAL conds = 0.0;
    char adumma[1];
    INTEGER idumma[1];
    INTEGER iwk = 0;
    INTEGER nnwork = 0;
    INTEGER isort = 0;
    char sort;
    INTEGER rsub = 0;
    INTEGER sdim = 0;
    INTEGER i = 0;
    INTEGER lwork = 0;
    REAL res[2];
    REAL tmp = 0.0;
    INTEGER knteig = 0;
    INTEGER ntest = 0;
    INTEGER nfail = 0;
    static const char *format_9992 = "(' Rdrves: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
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
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    path[0] = 'R';
    path[1] = 'E';
    path[2] = 'S';
    //
    //     Check for errors
    //
    ntestt = 0;
    ntestf = 0;
    info = 0;
    selopt = 0;
    //
    //     Important constants
    //
    badnn = false;
    nmax = 0;
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
    } else if (thresh < zero) {
        info = -6;
    } else if (nounit <= 0) {
        info = -7;
    } else if (lda < 1 || lda < nmax) {
        info = -9;
    } else if (ldvs < 1 || ldvs < nmax) {
        info = -17;
    } else if (5 * nmax + 2 * pow2(nmax) > nwork) {
        info = -20;
    }
    //
    if (info != 0) {
        Mxerbla("Rdrves", -info);
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
    ovfl = one / unfl;
    Rlabad(unfl, ovfl);
    ulp = Rlamch("Precision");
    ulpinv = one / ulp;
    rtulp = sqrt(ulp);
    rtulpi = one / rtulp;
    //
    //     Loop over sizes, types
    //
    nerrs = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        mtypes = maxtyp;
        if (nsizes == 1 && ntypes == maxtyp + 1) {
            mtypes++;
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_260;
            }
            //
            //           Save ISEED in case of an error.
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            //           Compute "A"
            //
            //           Control parameters:
            //
            //           KMAGN  KCONDS  KMODE        KTYPE
            //       =1  O(1)   1       clustered 1  zero
            //       =2  large  large   clustered 2  identity
            //       =3  small          exponential  Jordan
            //       =4                 arithmetic   diagonal, (w/ eigenvalues)
            //       =5                 random log   symmetric, w/ eigenvalues
            //       =6                 random       general, w/ eigenvalues
            //       =7                              random diagonal
            //       =8                              random symmetric
            //       =9                              random general
            //       =10                             random triangular
            //
            if (mtypes > maxtyp) {
                goto statement_90;
            }
            //
            itype = ktype[jtype - 1];
            imode = kmode[jtype - 1];
            //
            //           Compute norm
            //
            switch (kmagn[jtype - 1]) {
            case 1:
                goto statement_30;
            case 2:
                goto statement_40;
            case 3:
                goto statement_50;
            default:
                break;
            }
        //
        statement_30:
            anorm = one;
            goto statement_60;
        //
        statement_40:
            anorm = ovfl * ulp;
            goto statement_60;
        //
        statement_50:
            anorm = unfl * ulpinv;
            goto statement_60;
        //
        statement_60:
            //
            Rlaset("Full", lda, n, zero, zero, a, lda);
            iinfo = 0;
            cond = ulpinv;
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
                for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                    a[(jcol - 1) + (jcol - 1) * lda] = anorm;
                }
                //
            } else if (itype == 3) {
                //
                //              Jordan Block
                //
                for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                    a[(jcol - 1) + (jcol - 1) * lda] = anorm;
                    if (jcol > 1) {
                        a[(jcol - 1) + ((jcol - 1) - 1) * lda] = one;
                    }
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
            } else if (itype == 6) {
                //
                //              General, eigenvalues specified
                //
                if (kconds[jtype - 1] == 1) {
                    conds = one;
                } else if (kconds[jtype - 1] == 2) {
                    conds = rtulpi;
                } else {
                    conds = zero;
                }
                //
                adumma[1 - 1] = ' ';
                Rlatme(n, "S", iseed, work, imode, cond, one, adumma, "T", "T", "T", &work[(n + 1) - 1], 4, conds, n, n, anorm, a, lda, &work[(2 * n + 1) - 1], iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "S", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Symmetric, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "S", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              General, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                if (n >= 4) {
                    Rlaset("Full", 2, n, zero, zero, a, lda);
                    Rlaset("Full", n - 3, 1, zero, zero, &a[(3 - 1)], lda);
                    Rlaset("Full", n - 3, 2, zero, zero, &a[(3 - 1) + ((n - 1) - 1) * lda], lda);
                    Rlaset("Full", 1, n, zero, zero, &a[(n - 1)], lda);
                }
                //
            } else if (itype == 10) {
                //
                //              Triangular, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else {
                //
                iinfo = 1;
            }
            //
            if (iinfo != 0) {
                write(nounit, format_9992), "Generator", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                return;
            }
        //
        statement_90:
            //
            //           Test for minimal and generous workspace
            //
            for (iwk = 1; iwk <= 2; iwk = iwk + 1) {
                if (iwk == 1) {
                    nnwork = 3 * n;
                } else {
                    nnwork = 5 * n + 2 * pow2(n);
                }
                nnwork = max(nnwork, 1);
                //
                //              Initialize RESULT
                //
                for (j = 1; j <= 13; j = j + 1) {
                    result[j - 1] = -one;
                }
                //
                //              Test with and without sorting of eigenvalues
                //
                for (isort = 0; isort <= 1; isort = isort + 1) {
                    if (isort == 0) {
                        sort = 'N';
                        rsub = 0;
                    } else {
                        sort = 'S';
                        rsub = 6;
                    }
                    //
                    //                 Compute Schur form and Schur vectors, and test them
                    //
                    Rlacpy("F", n, n, a, lda, h, lda);
                    Rgees("V", &sort, Rslect, n, h, lda, sdim, wr, wi, vs, ldvs, work, nnwork, bwork, iinfo);
                    if (iinfo != 0 && iinfo != n + 2) {
                        result[(1 + rsub) - 1] = ulpinv;
                        write(nounit, format_9992), "Rgees1", iinfo, n, jtype, ioldsd;
                        info = abs(iinfo);
                        goto statement_220;
                    }
                    //
                    //                 Do Test (1) or Test (7)
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
                    //                 Do Tests (2) and (3) or Tests (8) and (9)
                    //
                    lwork = max((INTEGER)1, 2 * n * n);
                    Rhst01(n, 1, n, a, lda, h, lda, vs, ldvs, work, lwork, res);
                    result[(2 + rsub) - 1] = res[1 - 1];
                    result[(3 + rsub) - 1] = res[2 - 1];
                    //
                    //                 Do Test (4) or Test (10)
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
                            result[(4 + rsub) - 1] = max({result[(4 + rsub) - 1], abs(wi[i - 1] - tmp) / max(ulp * tmp, unfl)});
                            result[(4 + rsub) - 1] = max({result[(4 + rsub) - 1], abs(wi[(i + 1) - 1] + tmp) / max(ulp * tmp, unfl)});
                        } else if (i > 1) {
                            if (h[((i + 1) - 1) + (i - 1) * ldh] == zero && h[(i - 1) + ((i - 1) - 1) * ldh] == zero && wi[i - 1] != zero) {
                                result[(4 + rsub) - 1] = ulpinv;
                            }
                        }
                    }
                    //
                    //                 Do Test (5) or Test (11)
                    //
                    Rlacpy("F", n, n, a, lda, ht, lda);
                    Rgees("N", &sort, Rslect, n, ht, lda, sdim, wrt, wit, vs, ldvs, work, nnwork, bwork, iinfo);
                    if (iinfo != 0 && iinfo != n + 2) {
                        result[(5 + rsub) - 1] = ulpinv;
                        write(nounit, format_9992), "Rgees2", iinfo, n, jtype, ioldsd;
                        info = abs(iinfo);
                        goto statement_220;
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
                    //                 Do Test (6) or Test (12)
                    //
                    result[(6 + rsub) - 1] = zero;
                    for (i = 1; i <= n; i = i + 1) {
                        if (wr[i - 1] != wrt[i - 1] || wi[i - 1] != wit[i - 1]) {
                            result[(6 + rsub) - 1] = ulpinv;
                        }
                    }
                    //
                    //                 Do Test (13)
                    //
                    if (isort == 1) {
                        result[13 - 1] = zero;
                        knteig = 0;
                        for (i = 1; i <= n; i = i + 1) {
                            if (Rslect(wr[i - 1], wi[i - 1]) || Rslect(wr[i - 1], -wi[i - 1])) {
                                knteig++;
                            }
                            if (i < n) {
                                if ((Rslect(wr[(i + 1) - 1], wi[(i + 1) - 1]) || Rslect(wr[(i + 1) - 1], -wi[(i + 1) - 1])) && (!(Rslect(wr[i - 1], wi[i - 1]) || Rslect(wr[i - 1], -wi[i - 1]))) && iinfo != n + 2) {
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
            //              End of Loop -- Check for RESULT(j) > THRESH
            //
            statement_220:
                //
                ntest = 0;
                nfail = 0;
                for (j = 1; j <= 13; j = j + 1) {
                    if (result[j - 1] >= zero) {
                        ntest++;
                    }
                    if (result[j - 1] >= thresh) {
                        nfail++;
                    }
                }
                //
                if (nfail > 0) {
                    ntestf++;
                }
                if (ntestf == 1) {
                    write(nounit, "(/,1x,a3,' -- Real Schur Form Decomposition Driver',/,"
                                  "' Matrix types (see Rdrves for details): ')"),
                        path;
                    write(nounit, "(/,' Special Matrices:',/,'  1=Zero matrix.             ',"
                                  "'           ','  5=Diagonal: geometr. spaced entries.',/,"
                                  "'  2=Identity matrix.                    ','  6=Diagona',"
                                  "'l: clustered entries.',/,'  3=Transposed Jordan block.  ',"
                                  "'          ','  7=Diagonal: large, evenly spaced.',/,'  ',"
                                  "'4=Diagonal: evenly spaced entries.    ','  8=Diagonal: s',"
                                  "'mall, evenly spaced.')");
                    write(nounit, "(' Dense, Non-Symmetric Matrices:',/,'  9=Well-cond., ev',"
                                  "'enly spaced eigenvals.',' 14=Ill-cond., geomet. spaced e',"
                                  "'igenals.',/,' 10=Well-cond., geom. spaced eigenvals. ',"
                                  "' 15=Ill-conditioned, clustered e.vals.',/,' 11=Well-cond',"
                                  "'itioned, clustered e.vals. ',' 16=Ill-cond., random comp',"
                                  "'lex ',/,' 12=Well-cond., random complex ',6x,'   ',"
                                  "' 17=Ill-cond., large rand. complx ',/,' 13=Ill-condi',"
                                  "'tioned, evenly spaced.     ',' 18=Ill-cond., small rand.',"
                                  "' complx ')");
                    write(nounit, "(' 19=Matrix with random O(1) entries.    ',' 21=Matrix ',"
                                  "'with small random entries.',/,' 20=Matrix with large ran',"
                                  "'dom entries.   ',/)");
                    sprintnum_short(buf, thresh);
                    write(nounit, "(' Tests performed with test threshold =',a,/,"
                                  "' ( A denotes A on input and T denotes A on output)',/,/,"
                                  "' 1 = 0 if T in Schur form (no sort), ','  1/ulp otherwise',/,"
                                  "' 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)',/,"
                                  "' 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) ',/,"
                                  "' 4 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (no sort),',"
                                  "'  1/ulp otherwise',/,"
                                  "' 5 = 0 if T same no matter if VS computed (no sort),',"
                                  "'  1/ulp otherwise',/,"
                                  "' 6 = 0 if WR, WI same no matter if VS computed (no sort)',"
                                  "',  1/ulp otherwise')"),
                        buf;
                    write(nounit, "(' 7 = 0 if T in Schur form (sort), ','  1/ulp otherwise',/,"
                                  "' 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)',/,"
                                  "' 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) ',/,"
                                  "' 10 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (sort),',"
                                  "'  1/ulp otherwise',/,"
                                  "' 11 = 0 if T same no matter if VS computed (sort),',"
                                  "'  1/ulp otherwise',/,"
                                  "' 12 = 0 if WR, WI same no matter if VS computed (sort),',"
                                  "'  1/ulp otherwise',/,"
                                  "' 13 = 0 if sorting successful, 1/ulp otherwise',/)");
                    ntestf = 2;
                }
                //
                for (j = 1; j <= 13; j = j + 1) {
                    if (result[j - 1] >= thresh) {
                        sprintnum_short(buf, result[j - 1]);
                        write(nounit, "(' N=',i5,', IWK=',i2,', seed=',4(i4,','),' type ',i2,"
                                      "', test(',i2,')=',a)"),
                            n, iwk, ioldsd, jtype, j, buf;
                    }
                }
                //
                nerrs += nfail;
                ntestt += ntest;
                //
            }
        statement_260:;
        }
    }
    //
    //     Summary
    //
    Rlasum(path, nounit, nerrs, ntestt);
    //
    //     End of Rdrves
    //
}

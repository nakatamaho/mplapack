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

void Cdrgsx(INTEGER const nsize, INTEGER const ncmax, REAL const thresh, INTEGER const nin, INTEGER const nout, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *ai, COMPLEX *bi, COMPLEX *z, COMPLEX *q, COMPLEX *alpha, COMPLEX *beta, COMPLEX *c, INTEGER const ldc, REAL *s, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER const liwork, bool *bwork, INTEGER &info) {
    a([lda * star]);
    b([lda * star]);
    ai([lda * star]);
    bi([lda * star]);
    z([lda * star]);
    q([lda * star]);
    c([ldc * star]);
    common_read read(cmn);
    common_write write(cmn);
    INTEGER &m = cmn.m;
    INTEGER &n = cmn.n;
    INTEGER &mplusn = cmn.mplusn;
    INTEGER &k = cmn.k;
    bool &fs = cmn.fs;
    //
    COMPLEX x = 0.0;
    const REAL zero = 0.0;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER bdspac = 0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL ulpinv = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    const REAL ten = 1.0e+1;
    REAL thrsh2 = 0.0;
    INTEGER ntestt = 0;
    INTEGER nerrs = 0;
    INTEGER ifunc = 0;
    INTEGER prtype = 0;
    INTEGER qba = 0;
    INTEGER qbb = 0;
    REAL weight = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    char sense;
    INTEGER mm = 0;
    REAL pl[2];
    REAL difest[2];
    INTEGER linfo = 0;
    REAL result[10];
    REAL abnrm = 0.0;
    INTEGER ntest = 0;
    REAL temp1 = 0.0;
    INTEGER j = 0;
    bool ilabad = false;
    REAL temp2 = 0.0;
    INTEGER mn2 = 0;
    REAL diftru = 0.0;
    INTEGER i = 0;
    INTEGER nptknt = 0;
    REAL pltru = 0.0;
    static const char *format_9993 = "(/,' Tests performed:  (S is Schur, T is triangular, ','Q and Z are ',a,"
                                     "',',/,19x,' a is alpha, b is beta, and ',a,' means ',a,'.)',/,"
                                     "'  1 = | A - Q S Z',a,' | / ( |A| n ulp )      2 = | B - Q T Z',a,"
                                     "' | / ( |B| n ulp )',/,'  3 = | I - QQ',a,"
                                     "' | / ( n ulp )             4 = | I - ZZ',a,' | / ( n ulp )',/,"
                                     "'  5 = 1/ULP  if A is not in ','Schur form S',/,"
                                     "'  6 = difference between (alpha,beta)',' and diagonals of (S,T)',/,"
                                     "'  7 = 1/ULP  if SDIM is not the correct number of ',"
                                     "'selected eigenvalues',/,'  8 = 1/ULP  if DIFEST/DIFTRU > 10*THRESH or ',"
                                     "'DIFTRU/DIFEST > 10*THRESH',/,"
                                     "'  9 = 1/ULP  if DIFEST <> 0 or DIFTRU > ULP*norm(A,B) ',"
                                     "'when reordering fails',/,' 10 = 1/ULP  if PLEST/PLTRU > THRESH or ',"
                                     "'PLTRU/PLEST > THRESH',/,'    ( Test 10 is only for input examples )',/)";
    static const char *format_9996 = "(/,1x,a3,' -- Complex Expert Generalized Schur form',' problem driver')";
    static const char *format_9997 = "(' Cdrgsx: S not in Schur form at eigenvalue ',i6,'.',/,9x,'N=',i6,"
                                     "', JTYPE=',i6,')')";
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1(x) = abs(x.real()) + abs(x.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Check for errors
    //
    info = 0;
    if (nsize < 0) {
        info = -1;
    } else if (thresh < zero) {
        info = -2;
    } else if (nin <= 0) {
        info = -3;
    } else if (nout <= 0) {
        info = -4;
    } else if (lda < 1 || lda < nsize) {
        info = -6;
    } else if (ldc < 1 || ldc < nsize * nsize / 2) {
        info = -15;
    } else if (liwork < nsize + 2) {
        info = -21;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that point in the code,
    //       as well as the preferred amount for good performance.
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv.)
    //
    minwrk = 1;
    if (info == 0 && lwork >= 1) {
        minwrk = 3 * nsize * nsize / 2;
        //
        //        workspace for cggesx
        //
        maxwrk = nsize * (1 + iMlaenv(1, "Cgeqrf", " ", nsize, 1, nsize, 0));
        maxwrk = max({maxwrk, nsize * (1 + iMlaenv(1, "Cungqr", " ", nsize, 1, nsize, -1))});
        //
        //        workspace for Cgesvd
        //
        bdspac = 3 * nsize * nsize / 2;
        maxwrk = max({maxwrk, nsize * nsize * (1 + iMlaenv(1, "Cgebrd", " ", nsize * nsize / 2, nsize * nsize / 2, -1, -1))});
        maxwrk = max(maxwrk, bdspac);
        //
        maxwrk = max(maxwrk, minwrk);
        //
        work[1 - 1] = maxwrk;
    }
    //
    if (lwork < minwrk) {
        info = -18;
    }
    //
    if (info != 0) {
        Mxerbla("Cdrgsx", -info);
        return;
    }
    //
    //     Important constants
    //
    ulp = Rlamch("P");
    ulpinv = one / ulp;
    smlnum = Rlamch("S") / ulp;
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    thrsh2 = ten * thresh;
    ntestt = 0;
    nerrs = 0;
    //
    //     Go to the tests for read-in matrix pairs
    //
    ifunc = 0;
    if (nsize == 0) {
        goto statement_70;
    }
    //
    //     Test the built-in matrix pairs.
    //     Loop over different functions (IFUNC) of Cggesx, types (PRTYPE)
    //     of test matrices, different size (M+N)
    //
    prtype = 0;
    qba = 3;
    qbb = 4;
    weight = sqrt(ulp);
    //
    for (ifunc = 0; ifunc <= 3; ifunc = ifunc + 1) {
        for (prtype = 1; prtype <= 5; prtype = prtype + 1) {
            for (m = 1; m <= nsize - 1; m = m + 1) {
                for (n = 1; n <= nsize - m; n = n + 1) {
                    //
                    weight = one / weight;
                    mplusn = m + n;
                    //
                    //                 Generate test matrices
                    //
                    fs = true;
                    k = 0;
                    //
                    Claset("Full", mplusn, mplusn, czero, czero, ai, lda);
                    Claset("Full", mplusn, mplusn, czero, czero, bi, lda);
                    //
                    zlatm5(prtype, m, n, ai, lda, ai[((m + 1) - 1) + ((m + 1) - 1) * ldai], lda, ai[((m + 1) - 1) * ldai], lda, bi, lda, bi[((m + 1) - 1) + ((m + 1) - 1) * ldbi], lda, bi[((m + 1) - 1) * ldbi], lda, q, lda, z, lda, weight, qba, qbb);
                    //
                    //                 Compute the Schur factorization and swapping the
                    //                 m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
                    //                 Swapping is accomplished via the function Clctsx
                    //                 which is supplied below.
                    //
                    if (ifunc == 0) {
                        sense = "N";
                    } else if (ifunc == 1) {
                        sense = "E";
                    } else if (ifunc == 2) {
                        sense = "V";
                    } else if (ifunc == 3) {
                        sense = "B";
                    }
                    //
                    Clacpy("Full", mplusn, mplusn, ai, lda, a, lda);
                    Clacpy("Full", mplusn, mplusn, bi, lda, b, lda);
                    //
                    Cggesx("V", "V", "S", Clctsx, sense, mplusn, ai, lda, bi, lda, mm, alpha, beta, q, lda, z, lda, pl, difest, work, lwork, rwork, iwork, liwork, bwork, linfo);
                    //
                    if (linfo != 0 && linfo != mplusn + 2) {
                        result[1 - 1] = ulpinv;
                        write(nout, "(' Cdrgsx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
                                    "', JTYPE=',i6,')')"),
                            "Cggesx", linfo, mplusn, prtype;
                        info = linfo;
                        goto statement_30;
                    }
                    //
                    //                 Compute the norm(A, B)
                    //
                    Clacpy("Full", mplusn, mplusn, ai, lda, work, mplusn);
                    Clacpy("Full", mplusn, mplusn, bi, lda, &work[(mplusn * mplusn + 1) - 1], mplusn);
                    abnrm = Clange("Fro", mplusn, 2 * mplusn, work, mplusn, rwork);
                    //
                    //                 Do tests (1) to (4)
                    //
                    result[2 - 1] = zero;
                    Cget51(1, mplusn, a, lda, ai, lda, q, lda, z, lda, work, rwork, result[1 - 1]);
                    Cget51(1, mplusn, b, lda, bi, lda, q, lda, z, lda, work, rwork, result[2 - 1]);
                    Cget51(3, mplusn, b, lda, bi, lda, q, lda, q, lda, work, rwork, result[3 - 1]);
                    Cget51(3, mplusn, b, lda, bi, lda, z, lda, z, lda, work, rwork, result[4 - 1]);
                    ntest = 4;
                    //
                    //                 Do tests (5) and (6): check Schur form of A and
                    //                 compare eigenvalues with diagonals.
                    //
                    temp1 = zero;
                    result[5 - 1] = zero;
                    result[6 - 1] = zero;
                    //
                    for (j = 1; j <= mplusn; j = j + 1) {
                        ilabad = false;
                        temp2 = (abs1(alpha[j - 1] - ai[(j - 1) + (j - 1) * ldai]) / max({smlnum, abs1(alpha[j - 1]), abs1(ai[(j - 1) + (j - 1) * ldai])}) + abs1(beta[j - 1] - bi[(j - 1) + (j - 1) * ldbi]) / max({smlnum, abs1(beta[j - 1]), abs1(bi[(j - 1) + (j - 1) * ldbi])})) / ulp;
                        if (j < mplusn) {
                            if (ai[((j + 1) - 1) + (j - 1) * ldai] != zero) {
                                ilabad = true;
                                result[5 - 1] = ulpinv;
                            }
                        }
                        if (j > 1) {
                            if (ai[(j - 1) + ((j - 1) - 1) * ldai] != zero) {
                                ilabad = true;
                                result[5 - 1] = ulpinv;
                            }
                        }
                        temp1 = max(temp1, temp2);
                        if (ilabad) {
                            write(nout, format_9997), j, mplusn, prtype;
                        }
                    }
                    result[6 - 1] = temp1;
                    ntest += 2;
                    //
                    //                 Test (7) (if sorting worked)
                    //
                    result[7 - 1] = zero;
                    if (linfo == mplusn + 3) {
                        result[7 - 1] = ulpinv;
                    } else if (mm != n) {
                        result[7 - 1] = ulpinv;
                    }
                    ntest++;
                    //
                    //                 Test (8): compare the estimated value DIF and its
                    //                 value. first, compute the exact DIF.
                    //
                    result[8 - 1] = zero;
                    mn2 = mm * (mplusn - mm) * 2;
                    if (ifunc >= 2 && mn2 <= ncmax * ncmax) {
                        //
                        //                    Note: for either following two cases, there are
                        //                    almost same number of test cases fail the test.
                        //
                        zlakf2(mm, mplusn - mm, ai, lda, ai[((mm + 1) - 1) + ((mm + 1) - 1) * ldai], bi, bi[((mm + 1) - 1) + ((mm + 1) - 1) * ldbi], c, ldc);
                        //
                        Cgesvd("N", "N", mn2, mn2, c, ldc, s, work, 1, &work[2 - 1], 1, &work[3 - 1], lwork - 2, rwork, info);
                        diftru = s[mn2 - 1];
                        //
                        if (difest[2 - 1] == zero) {
                            if (diftru > abnrm * ulp) {
                                result[8 - 1] = ulpinv;
                            }
                        } else if (diftru == zero) {
                            if (difest[2 - 1] > abnrm * ulp) {
                                result[8 - 1] = ulpinv;
                            }
                        } else if ((diftru > thrsh2 * difest[2 - 1]) || (diftru * thrsh2 < difest[2 - 1])) {
                            result[8 - 1] = max(diftru / difest[2 - 1], difest[2 - 1] / diftru);
                        }
                        ntest++;
                    }
                    //
                    //                 Test (9)
                    //
                    result[9 - 1] = zero;
                    if (linfo == (mplusn + 2)) {
                        if (diftru > abnrm * ulp) {
                            result[9 - 1] = ulpinv;
                        }
                        if ((ifunc > 1) && (difest[2 - 1] != zero)) {
                            result[9 - 1] = ulpinv;
                        }
                        if ((ifunc == 1) && (pl[1 - 1] != zero)) {
                            result[9 - 1] = ulpinv;
                        }
                        ntest++;
                    }
                    //
                    ntestt += ntest;
                    //
                    //                 Print out tests which fail.
                    //
                    for (j = 1; j <= 9; j = j + 1) {
                        if (result[j - 1] >= thresh) {
                            //
                            //                       If this is the first test to fail,
                            //                       prINTEGER a header to the data file.
                            //
                            if (nerrs == 0) {
                                write(nout, format_9996), "ZGX";
                                //
                                //                          Matrix types
                                //
                                write(nout, "(' Matrix types: ',/,"
                                            "'  1:  A is a block diagonal matrix of Jordan blocks ',"
                                            "'and B is the identity ',/,'      matrix, ',/,"
                                            "'  2:  A and B are upper triangular matrices, ',/,"
                                            "'  3:  A and B are as type 2, but each second diagonal ',"
                                            "'block in A_11 and ',/,"
                                            "'      each third diaongal block in A_22 are 2x2 blocks,',"
                                            "/,'  4:  A and B are block diagonal matrices, ',/,"
                                            "'  5:  (A,B) has potentially close or common ',"
                                            "'eigenvalues.',/)");
                                //
                                //                          Tests performed
                                //
                                {
                                    write_loop wloop(cmn, nout, format_9993);
                                    wloop, "unitary", "'", "transpose";
                                    for (i = 1; i <= 4; i = i + 1) {
                                        wloop, "'";
                                    }
                                }
                                //
                            }
                            nerrs++;
                            if (result[j - 1] < 10000.0) {
                                write(nout, "(' Matrix order=',i2,', type=',i2,', a=',d10.3,"
                                            "', order(A_11)=',i2,', result ',i2,' is ',0p,f8.2)"),
                                    mplusn, prtype, weight, m, j, result(j);
                            } else {
                                write(nout, "(' Matrix order=',i2,', type=',i2,', a=',d10.3,"
                                            "', order(A_11)=',i2,', result ',i2,' is ',0p,d10.3)"),
                                    mplusn, prtype, weight, m, j, result(j);
                            }
                        }
                    }
                //
                statement_30:;
                }
            }
        }
    }
    //
    goto statement_150;
//
statement_70:
    //
    //     Read in data from file to check accuracy of condition estimation
    //     Read input data until N=0
    //
    nptknt = 0;
//
statement_80:
    try {
        read(nin, star), mplusn;
    } catch (read_end const) {
        goto statement_140;
    }
    if (mplusn == 0) {
        goto statement_140;
    }
    try {
        read(nin, star), n;
    } catch (read_end const) {
        goto statement_140;
    }
    for (i = 1; i <= mplusn; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= mplusn; j = j + 1) {
                rloop, ai(i, j);
            }
        }
    }
    for (i = 1; i <= mplusn; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= mplusn; j = j + 1) {
                rloop, bi(i, j);
            }
        }
    }
    read(nin, star), pltru, diftru;
    //
    nptknt++;
    fs = true;
    k = 0;
    m = mplusn - n;
    //
    Clacpy("Full", mplusn, mplusn, ai, lda, a, lda);
    Clacpy("Full", mplusn, mplusn, bi, lda, b, lda);
    //
    //     Compute the Schur factorization while swapping the
    //     m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
    //
    Cggesx("V", "V", "S", Clctsx, "B", mplusn, ai, lda, bi, lda, mm, alpha, beta, q, lda, z, lda, pl, difest, work, lwork, rwork, iwork, liwork, bwork, linfo);
    //
    if (linfo != 0 && linfo != mplusn + 2) {
        result[1 - 1] = ulpinv;
        write(nout, "(' Cdrgsx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
                    "', Input Example #',i2,')')"),
            "Cggesx", linfo, mplusn, nptknt;
        goto statement_130;
    }
    //
    //     Compute the norm(A, B)
    //        (should this be norm of (A,B) or (AI,BI)?)
    //
    Clacpy("Full", mplusn, mplusn, ai, lda, work, mplusn);
    Clacpy("Full", mplusn, mplusn, bi, lda, &work[(mplusn * mplusn + 1) - 1], mplusn);
    abnrm = Clange("Fro", mplusn, 2 * mplusn, work, mplusn, rwork);
    //
    //     Do tests (1) to (4)
    //
    Cget51(1, mplusn, a, lda, ai, lda, q, lda, z, lda, work, rwork, result[1 - 1]);
    Cget51(1, mplusn, b, lda, bi, lda, q, lda, z, lda, work, rwork, result[2 - 1]);
    Cget51(3, mplusn, b, lda, bi, lda, q, lda, q, lda, work, rwork, result[3 - 1]);
    Cget51(3, mplusn, b, lda, bi, lda, z, lda, z, lda, work, rwork, result[4 - 1]);
    //
    //     Do tests (5) and (6): check Schur form of A and compare
    //     eigenvalues with diagonals.
    //
    ntest = 6;
    temp1 = zero;
    result[5 - 1] = zero;
    result[6 - 1] = zero;
    //
    for (j = 1; j <= mplusn; j = j + 1) {
        ilabad = false;
        temp2 = (abs1(alpha[j - 1] - ai[(j - 1) + (j - 1) * ldai]) / max({smlnum, abs1(alpha[j - 1]), abs1(ai[(j - 1) + (j - 1) * ldai])}) + abs1(beta[j - 1] - bi[(j - 1) + (j - 1) * ldbi]) / max({smlnum, abs1(beta[j - 1]), abs1(bi[(j - 1) + (j - 1) * ldbi])})) / ulp;
        if (j < mplusn) {
            if (ai[((j + 1) - 1) + (j - 1) * ldai] != zero) {
                ilabad = true;
                result[5 - 1] = ulpinv;
            }
        }
        if (j > 1) {
            if (ai[(j - 1) + ((j - 1) - 1) * ldai] != zero) {
                ilabad = true;
                result[5 - 1] = ulpinv;
            }
        }
        temp1 = max(temp1, temp2);
        if (ilabad) {
            write(nout, format_9997), j, mplusn, nptknt;
        }
    }
    result[6 - 1] = temp1;
    //
    //     Test (7) (if sorting worked)  <--------- need to be checked.
    //
    ntest = 7;
    result[7 - 1] = zero;
    if (linfo == mplusn + 3) {
        result[7 - 1] = ulpinv;
    }
    //
    //     Test (8): compare the estimated value of DIF and its true value.
    //
    ntest = 8;
    result[8 - 1] = zero;
    if (difest[2 - 1] == zero) {
        if (diftru > abnrm * ulp) {
            result[8 - 1] = ulpinv;
        }
    } else if (diftru == zero) {
        if (difest[2 - 1] > abnrm * ulp) {
            result[8 - 1] = ulpinv;
        }
    } else if ((diftru > thrsh2 * difest[2 - 1]) || (diftru * thrsh2 < difest[2 - 1])) {
        result[8 - 1] = max(diftru / difest[2 - 1], difest[2 - 1] / diftru);
    }
    //
    //     Test (9)
    //
    ntest = 9;
    result[9 - 1] = zero;
    if (linfo == (mplusn + 2)) {
        if (diftru > abnrm * ulp) {
            result[9 - 1] = ulpinv;
        }
        if ((ifunc > 1) && (difest[2 - 1] != zero)) {
            result[9 - 1] = ulpinv;
        }
        if ((ifunc == 1) && (pl[1 - 1] != zero)) {
            result[9 - 1] = ulpinv;
        }
    }
    //
    //     Test (10): compare the estimated value of PL and it true value.
    //
    ntest = 10;
    result[10 - 1] = zero;
    if (pl[1 - 1] == zero) {
        if (pltru > abnrm * ulp) {
            result[10 - 1] = ulpinv;
        }
    } else if (pltru == zero) {
        if (pl[1 - 1] > abnrm * ulp) {
            result[10 - 1] = ulpinv;
        }
    } else if ((pltru > thresh * pl[1 - 1]) || (pltru * thresh < pl[1 - 1])) {
        result[10 - 1] = ulpinv;
    }
    //
    ntestt += ntest;
    //
    //     Print out tests which fail.
    //
    for (j = 1; j <= ntest; j = j + 1) {
        if (result[j - 1] >= thresh) {
            //
            //           If this is the first test to fail,
            //           prINTEGER a header to the data file.
            //
            if (nerrs == 0) {
                write(nout, format_9996), "ZGX";
                //
                //              Matrix types
                //
                write(nout, "('Input Example')");
                //
                //              Tests performed
                //
                {
                    write_loop wloop(cmn, nout, format_9993);
                    wloop, "unitary", "'", "transpose";
                    for (i = 1; i <= 4; i = i + 1) {
                        wloop, "'";
                    }
                }
                //
            }
            nerrs++;
            if (result[j - 1] < 10000.0) {
                write(nout, "(' Input example #',i2,', matrix order=',i4,',',' result ',i2,"
                            "' is',0p,f8.2)"),
                    nptknt, mplusn, j, result(j);
            } else {
                write(nout, "(' Input example #',i2,', matrix order=',i4,',',' result ',i2,"
                            "' is',1p,d10.3)"),
                    nptknt, mplusn, j, result(j);
            }
        }
        //
    }
//
statement_130:
    goto statement_80;
statement_140:
//
statement_150:
    //
    //     Summary
    //
    Alasvm("ZGX", nout, nerrs, ntestt, 0);
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Cdrgsx
    //
}

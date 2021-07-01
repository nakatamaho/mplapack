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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;
using std::regex;
using std::regex_replace;

void Rdrgsx(INTEGER const nsize, INTEGER const ncmax, REAL const thresh, INTEGER const nin, INTEGER const nout, REAL *a, INTEGER const lda, REAL *b, REAL *ai, REAL *bi, REAL *z, REAL *q, REAL *alphar, REAL *alphai, REAL *beta, REAL *c, INTEGER const ldc, REAL *s, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, bool *bwork, INTEGER &info) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    INTEGER m;
    INTEGER n;
    INTEGER mplusn;
    INTEGER k;
    bool fs;
    INTEGER ldai = lda;
    INTEGER ldbi = lda;
    INTEGER ldq = lda;
    char buf0[1024];
    char buf1[1024];
    double dtmp, dtmp1, dtmp2;
    //
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
    INTEGER i1 = 0;
    INTEGER iinfo = 0;
    INTEGER mn2 = 0;
    REAL diftru = 0.0;
    INTEGER i = 0;
    INTEGER nptknt = 0;
    REAL pltru = 0.0;
    string str;
    istringstream iss;

    static const char *format_9992 = "(/,' Tests performed:  (S is Schur, T is triangular, ','Q and Z are ',a,"
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
    static const char *format_9995 = "(/,1x,a3,' -- Real Expert Generalized Schur form',' problem driver')";
    static const char *format_9996 = "(' Rdrgsx: S not in Schur form at eigenvalue ',i6,'.',/,9x,'N=',i6,"
                                     "', JTYPE=',i6,')')";
    static const char *format_9997 = "(' Rdrgsx: Rget53 returned INFO=',i1,' for eigenvalue ',i6,'.',/,9x,'N=',"
                                     "i6,', JTYPE=',i6,')')";
    //
    //     Check for errors
    //
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
        info = -17;
    } else if (liwork < nsize + 6) {
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
        minwrk = max((INTEGER)10 * (nsize + 1), 5 * nsize * nsize / 2);
        //
        //        workspace for sggesx
        //
        maxwrk = 9 * (nsize + 1) + nsize * iMlaenv(1, "Rgeqrf", " ", nsize, 1, nsize, 0);
        maxwrk = max({maxwrk, 9 * (nsize + 1) + nsize * iMlaenv(1, "Rorgqr", " ", nsize, 1, nsize, -1)});
        //
        //        workspace for Rgesvd
        //
        bdspac = 5 * nsize * nsize / 2;
        maxwrk = max({maxwrk, 3 * nsize * nsize / 2 + nsize * nsize * iMlaenv(1, "Rgebrd", " ", nsize * nsize / 2, nsize * nsize / 2, -1, -1)});
        maxwrk = max(maxwrk, bdspac);
        //
        maxwrk = max(maxwrk, minwrk);
        //
        work[1 - 1] = maxwrk;
    }
    //
    if (lwork < minwrk) {
        info = -19;
    }
    //
    if (info != 0) {
        Mxerbla("Rdrgsx", -info);
        return;
    }
    //
    //     Important constants
    //
    ulp = Rlamch("P");
    ulpinv = one / ulp;
    smlnum = Rlamch("S") / ulp;
    bignum = one / smlnum;
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
    //     Loop over different functions (IFUNC) of Rggesx, types (PRTYPE)
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
                    Rlaset("Full", mplusn, mplusn, zero, zero, ai, lda);
                    Rlaset("Full", mplusn, mplusn, zero, zero, bi, lda);
                    //
                    Rlatm5(prtype, m, n, ai, lda, &ai[((m + 1) - 1) + ((m + 1) - 1) * ldai], lda, &ai[((m + 1) - 1) * ldai], lda, bi, lda, &bi[((m + 1) - 1) + ((m + 1) - 1) * ldbi], lda, &bi[((m + 1) - 1) * ldbi], lda, q, lda, z, lda, weight, qba, qbb);
                    //
                    //                 Compute the Schur factorization and swapping the
                    //                 m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
                    //                 Swapping is accomplished via the function Rlctsx
                    //                 which is supplied below.
                    //
                    if (ifunc == 0) {
                        sense = 'N';
                    } else if (ifunc == 1) {
                        sense = 'E';
                    } else if (ifunc == 2) {
                        sense = 'V';
                    } else if (ifunc == 3) {
                        sense = 'B';
                    }
                    //
                    Rlacpy("Full", mplusn, mplusn, ai, lda, a, lda);
                    Rlacpy("Full", mplusn, mplusn, bi, lda, b, lda);
                    //
                    Rggesx("V", "V", "S", Rlctsx, &sense, mplusn, ai, lda, bi, lda, mm, alphar, alphai, beta, q, lda, z, lda, pl, difest, work, lwork, iwork, liwork, bwork, linfo);
                    //
                    if (linfo != 0 && linfo != mplusn + 2) {
                        result[1 - 1] = ulpinv;
                        write(nout, "(' Rdrgsx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
                                    "', JTYPE=',i6,')')"),
                            "Rggesx", linfo, mplusn, prtype;
                        info = linfo;
                        goto statement_30;
                    }
                    //
                    //                 Compute the norm(A, B)
                    //
                    Rlacpy("Full", mplusn, mplusn, ai, lda, work, mplusn);
                    Rlacpy("Full", mplusn, mplusn, bi, lda, &work[(mplusn * mplusn + 1) - 1], mplusn);
                    abnrm = Rlange("Fro", mplusn, 2 * mplusn, work, mplusn, work);
                    //
                    //                 Do tests (1) to (4)
                    //
                    Rget51(1, mplusn, a, lda, ai, lda, q, lda, z, lda, work, result[1 - 1]);
                    Rget51(1, mplusn, b, lda, bi, lda, q, lda, z, lda, work, result[2 - 1]);
                    Rget51(3, mplusn, b, lda, bi, lda, q, lda, q, lda, work, result[3 - 1]);
                    Rget51(3, mplusn, b, lda, bi, lda, z, lda, z, lda, work, result[4 - 1]);
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
                        if (alphai[j - 1] == zero) {
                            temp2 = (abs(alphar[j - 1] - ai[(j - 1) + (j - 1) * ldai]) / max({smlnum, abs(alphar[j - 1]), abs(ai[(j - 1) + (j - 1) * ldai])}) + abs(beta[j - 1] - bi[(j - 1) + (j - 1) * ldbi]) / max({smlnum, abs(beta[j - 1]), abs(bi[(j - 1) + (j - 1) * ldbi])})) / ulp;
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
                        } else {
                            if (alphai[j - 1] > zero) {
                                i1 = j;
                            } else {
                                i1 = j - 1;
                            }
                            if (i1 <= 0 || i1 >= mplusn) {
                                ilabad = true;
                            } else if (i1 < mplusn - 1) {
                                if (ai[((i1 + 2) - 1) + ((i1 + 1) - 1) * ldai] != zero) {
                                    ilabad = true;
                                    result[5 - 1] = ulpinv;
                                }
                            } else if (i1 > 1) {
                                if (ai[(i1 - 1) + ((i1 - 1) - 1) * ldai] != zero) {
                                    ilabad = true;
                                    result[5 - 1] = ulpinv;
                                }
                            }
                            if (!ilabad) {
                                Rget53(&ai[(i1 - 1) + (i1 - 1) * ldai], lda, &bi[(i1 - 1) + (i1 - 1) * ldbi], lda, beta[j - 1], alphar[j - 1], alphai[j - 1], temp2, iinfo);
                                if (iinfo >= 3) {
                                    write(nout, format_9997), iinfo, j, mplusn, prtype;
                                    info = abs(iinfo);
                                }
                            } else {
                                temp2 = ulpinv;
                            }
                        }
                        temp1 = max(temp1, temp2);
                        if (ilabad) {
                            write(nout, format_9996), j, mplusn, prtype;
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
                        //                    Note: for either following two causes, there are
                        //                    almost same number of test cases fail the test.
                        //
                        Rlakf2(mm, mplusn - mm, ai, lda, &ai[((mm + 1) - 1) + ((mm + 1) - 1) * ldai], bi, &bi[((mm + 1) - 1) + ((mm + 1) - 1) * ldbi], c, ldc);
                        //
                        Rgesvd("N", "N", mn2, mn2, c, ldc, s, work, 1, &work[2 - 1], 1, &work[3 - 1], lwork - 2, info);
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
                            //                       print a header to the data file.
                            //
                            if (nerrs == 0) {
                                write(nout, format_9995), "DGX";
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
                                    write_loop wloop(cmn, nout, format_9992);
                                    wloop, "orthogonal", "'", "transpose";
                                    for (i = 1; i <= 4; i = i + 1) {
                                        wloop, "'";
                                    }
                                }
                                //
                            }
                            nerrs++;
                            if (result[j - 1] < 10000.0) {
                                sprintnum_short(buf0, weight);
                                sprintnum_short(buf1, result[j - 1]);
                                write(nout, "(' Matrix order=',i2,', type=',i2,', a=',a,"
                                            "', order(A_11)=',i2,', result ',i2,' is ',0p,a)"),
                                    mplusn, prtype, buf0, m, j, buf1;
                            } else {
                                sprintnum_short(buf0, weight);
                                sprintnum_short(buf1, result[j - 1]);
                                write(nout, "(' Matrix order=',i2,', type=',i2,', a=',a,"
                                            "', order(A_11)=',i2,', result ',i2,' is ',0p,a)"),
                                    mplusn, prtype, buf0, m, j, buf1;
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

    while (1) {
        getline(cin, str);
        iss.clear();
        iss.str(str);
        iss >> mplusn;
        if (mplusn == 0)
            break;
        getline(cin, str);
        iss.clear();
        iss.str(str);
        iss >> n;
        for (i = 1; i <= mplusn; i = i + 1) {
            getline(cin, str);
            string _r = regex_replace(str, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            for (j = 1; j <= mplusn; j = j + 1) {
                iss >> dtmp;
                ai[(i - 1) + (j - 1) * ldai] = dtmp;
            }
        }
        for (i = 1; i <= mplusn; i = i + 1) {
            getline(cin, str);
            string _r = regex_replace(str, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            for (j = 1; j <= mplusn; j = j + 1) {
                iss >> dtmp;
                bi[(i - 1) + (j - 1) * ldbi] = dtmp;
            }
        }
        getline(cin, str);
        string _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        iss.clear();
        iss.str(str);
        iss >> dtmp;
        pltru = dtmp;
        iss >> dtmp;
        diftru = dtmp;
        //
        nptknt++;
        fs = true;
        k = 0;
        m = mplusn - n;
        //
        Rlacpy("Full", mplusn, mplusn, ai, lda, a, lda);
        Rlacpy("Full", mplusn, mplusn, bi, lda, b, lda);
        //
        //     Compute the Schur factorization while swapping the
        //     m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
        //
        Rggesx("V", "V", "S", Rlctsx, "B", mplusn, ai, lda, bi, lda, mm, alphar, alphai, beta, q, lda, z, lda, pl, difest, work, lwork, iwork, liwork, bwork, linfo);
        //
        if (linfo != 0 && linfo != mplusn + 2) {
            result[1 - 1] = ulpinv;
            write(nout, "(' Rdrgsx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
                        "', Input Example #',i2,')')"),
                "Rggesx", linfo, mplusn, nptknt;
            continue;
        }
        //
        //     Compute the norm(A, B)
        //        (should this be norm of (A,B) or (AI,BI)?)
        //
        Rlacpy("Full", mplusn, mplusn, ai, lda, work, mplusn);
        Rlacpy("Full", mplusn, mplusn, bi, lda, &work[(mplusn * mplusn + 1) - 1], mplusn);
        abnrm = Rlange("Fro", mplusn, 2 * mplusn, work, mplusn, work);
        //
        //     Do tests (1) to (4)
        //
        Rget51(1, mplusn, a, lda, ai, lda, q, lda, z, lda, work, result[1 - 1]);
        Rget51(1, mplusn, b, lda, bi, lda, q, lda, z, lda, work, result[2 - 1]);
        Rget51(3, mplusn, b, lda, bi, lda, q, lda, q, lda, work, result[3 - 1]);
        Rget51(3, mplusn, b, lda, bi, lda, z, lda, z, lda, work, result[4 - 1]);
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
            if (alphai[j - 1] == zero) {
                temp2 = (abs(alphar[j - 1] - ai[(j - 1) + (j - 1) * ldai]) / max({smlnum, abs(alphar[j - 1]), abs(ai[(j - 1) + (j - 1) * ldai])}) + abs(beta[j - 1] - bi[(j - 1) + (j - 1) * ldbi]) / max({smlnum, abs(beta[j - 1]), abs(bi[(j - 1) + (j - 1) * ldbi])})) / ulp;
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
            } else {
                if (alphai[j - 1] > zero) {
                    i1 = j;
                } else {
                    i1 = j - 1;
                }
                if (i1 <= 0 || i1 >= mplusn) {
                    ilabad = true;
                } else if (i1 < mplusn - 1) {
                    if (ai[((i1 + 2) - 1) + ((i1 + 1) - 1) * ldai] != zero) {
                        ilabad = true;
                        result[5 - 1] = ulpinv;
                    }
                } else if (i1 > 1) {
                    if (ai[(i1 - 1) + ((i1 - 1) - 1) * ldai] != zero) {
                        ilabad = true;
                        result[5 - 1] = ulpinv;
                    }
                }
                if (!ilabad) {
                    Rget53(&ai[(i1 - 1) + (i1 - 1) * ldai], lda, &bi[(i1 - 1) + (i1 - 1) * ldbi], lda, beta[j - 1], alphar[j - 1], alphai[j - 1], temp2, iinfo);
                    if (iinfo >= 3) {
                        write(nout, format_9997), iinfo, j, mplusn, nptknt;
                        info = abs(iinfo);
                    }
                } else {
                    temp2 = ulpinv;
                }
            }
            temp1 = max(temp1, temp2);
            if (ilabad) {
                write(nout, format_9996), j, mplusn, nptknt;
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
                //           print a header to the data file.
                //
                if (nerrs == 0) {
                    write(nout, format_9995), "DGX";
                    //
                    //              Matrix types
                    //
                    write(nout, "('Input Example')");
                    //
                    //              Tests performed
                    //
                    {
                        write_loop wloop(cmn, nout, format_9992);
                        wloop, "orthogonal", "'", "transpose";
                        for (i = 1; i <= 4; i = i + 1) {
                            wloop, "'";
                        }
                    }
                    //
                }
                nerrs++;
                if (result[j - 1] < 10000.0) {
                    sprintnum_short(buf0, result[j - 1]);
                    write(nout, "(' Input example #',i2,', matrix order=',i4,',',' result ',i2,"
                                "' is',0p,a)"),
                        nptknt, mplusn, j, buf0;
                } else {
                    sprintnum_short(buf0, result[j - 1]);
                    write(nout, "(' Input example #',i2,', matrix order=',i4,',',' result ',i2,"
                                "' is',1p,a)"),
                        nptknt, mplusn, j, buf0;
                }
            }
            //
        }
        //
    }
statement_150:
    //
    //     Summary
    //
    Alasvm("DGX", nout, nerrs, ntestt, 0);
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Rdrgsx
    //
}

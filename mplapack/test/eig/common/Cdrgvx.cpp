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

void Cdrgvx(INTEGER const nsize, REAL const thresh, INTEGER const nin, INTEGER const nout, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *ai, COMPLEX *bi, COMPLEX *alpha, COMPLEX *beta, COMPLEX *vl, COMPLEX *vr, INTEGER const ilo, INTEGER const ihi, REAL *lscale, REAL *rscale, REAL *s, REAL *dtru, REAL *dif, REAL *diftru, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER const liwork, REAL *result, bool *bwork, INTEGER &info) {
    INTEGER ldb = lda;
    INTEGER ldai = lda;
    INTEGER ldbi = lda;
    INTEGER ldvl = lda;
    INTEGER ldvr = lda;
    char buf[1024];
    common cmn;
    common_write write(cmn);
    INTEGER nmax = 0;
    const REAL zero = 0.0;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER n = 0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL ulpinv = 0.0;
    const REAL ten = 1.0e+1;
    REAL thrsh2 = 0.0;
    INTEGER nerrs = 0;
    INTEGER nptknt = 0;
    INTEGER ntestt = 0;
    const REAL tnth = 1.0e-1;
    COMPLEX weight[5];
    const REAL half = 0.5e+0;
    INTEGER iptype = 0;
    INTEGER iwa = 0;
    INTEGER iwb = 0;
    INTEGER iwx = 0;
    INTEGER iwy = 0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER linfo = 0;
    REAL abnorm = 0.0;
    INTEGER i = 0;
    REAL ratio1 = 0.0;
    REAL ratio2 = 0.0;
    INTEGER j = 0;
    static const char *format_9986 = "(' Cdrgvx: ',a,' Eigenvectors from ',a,' incorrectly ','normalized.',/,"
                                     "' Bits of error=',0p,a,',',9x,'N=',i6,', Input Example #',i2,')')";
    static const char *format_9992 = "(/,' Tests performed:  ',/,4x,"
                                     "' a is alpha, b is beta, l is a left eigenvector, ',/,4x,"
                                     "' r is a right eigenvector and ',a,' means ',a,'.',/,"
                                     "' 1 = max | ( b A - a B )',a,' l | / const.',/,"
                                     "' 2 = max | ( b A - a B ) r | / const.',/,"
                                     "' 3 = max ( Sest/Stru, Stru/Sest ) ',' over all eigenvalues',/,"
                                     "' 4 = max( DIFest/DIFtru, DIFtru/DIFest ) ',"
                                     "' over the 1st and 5th eigenvectors',/)";
    static const char *format_9997 = "(/,1x,a3,' -- Complex Expert Eigenvalue/vector',' problem driver')";
    static const char *format_9998 = "(' Cdrgvx: ',a,' Eigenvectors from ',a,' incorrectly ','normalized.',/,"
                                     "' Bits of error=',0p,a,',',9x,'N=',i6,', JTYPE=',i6,', IWA=',i5,"
                                     "', IWB=',i5,', IWX=',i5,', IWY=',i5)";
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
    //     .. Executable Statements ..
    //
    //     Check for errors
    //
    info = 0;
    //
    nmax = 5;
    //
    if (nsize < 0) {
        info = -1;
    } else if (thresh < zero) {
        info = -2;
    } else if (nin <= 0) {
        info = -3;
    } else if (nout <= 0) {
        info = -4;
    } else if (lda < 1 || lda < nmax) {
        info = -6;
    } else if (liwork < nmax + 2) {
        info = -26;
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
        minwrk = 2 * nmax * (nmax + 1);
        maxwrk = nmax * (1 + iMlaenv(1, "Cgeqrf", " ", nmax, 1, nmax, 0));
        maxwrk = max(maxwrk, 2 * nmax * (nmax + 1));
        work[1 - 1] = maxwrk;
    }
    //
    if (lwork < minwrk) {
        info = -23;
    }
    //
    if (info != 0) {
        Mxerbla("Cdrgvx", -info);
        return;
    }
    //
    n = 5;
    ulp = Rlamch("P");
    ulpinv = one / ulp;
    thrsh2 = ten * thresh;
    nerrs = 0;
    nptknt = 0;
    ntestt = 0;
    //
    if (nsize == 0) {
        goto statement_90;
    }
    //
    //     Parameters used for generating test matrices.
    //
    weight[1 - 1] = COMPLEX(tnth, zero);
    weight[2 - 1] = COMPLEX(half, zero);
    weight[3 - 1] = one;
    weight[4 - 1] = one / weight[2 - 1];
    weight[5 - 1] = one / weight[1 - 1];
    //
    for (iptype = 1; iptype <= 2; iptype = iptype + 1) {
        for (iwa = 1; iwa <= 5; iwa = iwa + 1) {
            for (iwb = 1; iwb <= 5; iwb = iwb + 1) {
                for (iwx = 1; iwx <= 5; iwx = iwx + 1) {
                    for (iwy = 1; iwy <= 5; iwy = iwy + 1) {
                        //
                        //                    generated a pair of test matrix
                        //
                        Clatm6(iptype, 5, a, lda, b, vr, lda, vl, lda, weight[iwa - 1], weight[iwb - 1], weight[iwx - 1], weight[iwy - 1], dtru, diftru);
                        //
                        //                    Compute eigenvalues/eigenvectors of (A, B).
                        //                    Compute eigenvalue/eigenvector condition numbers
                        //                    using computed eigenvectors.
                        //
                        Clacpy("F", n, n, a, lda, ai, lda);
                        Clacpy("F", n, n, b, lda, bi, lda);
                        //
                        Cggevx("N", "V", "V", "B", n, ai, lda, bi, lda, alpha, beta, vl, lda, vr, lda, ilo, ihi, lscale, rscale, anorm, bnorm, s, dif, work, lwork, rwork, iwork, bwork, linfo);
                        if (linfo != 0) {
                            write(nout, "(' Cdrgvx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
                                        "', JTYPE=',i6,')')"),
                                "Cggevx", linfo, n, iptype, iwa, iwb, iwx, iwy;
                            goto statement_30;
                        }
                        //
                        //                    Compute the norm(A, B)
                        //
                        Clacpy("Full", n, n, ai, lda, work, n);
                        Clacpy("Full", n, n, bi, lda, &work[(n * n + 1) - 1], n);
                        abnorm = Clange("Fro", n, 2 * n, work, n, rwork);
                        //
                        //                    Tests (1) and (2)
                        //
                        result[1 - 1] = zero;
                        Cget52(true, n, a, lda, b, lda, vl, lda, alpha, beta, work, rwork, &result[1 - 1]);
                        if (result[2 - 1] > thresh) {
                            sprintnum_short(buf, result[2 - 1]);
                            write(nout, format_9998), "Left", "Cggevx", buf, n, iptype, iwa, iwb, iwx, iwy;
                        }
                        //
                        result[2 - 1] = zero;
                        Cget52(false, n, a, lda, b, lda, vr, lda, alpha, beta, work, rwork, &result[2 - 1]);
                        if (result[3 - 1] > thresh) {
                            sprintnum_short(buf, result[3 - 1]);
                            write(nout, format_9998), "Right", "Cggevx", buf, n, iptype, iwa, iwb, iwx, iwy;
                        }
                        //
                        //                    Test (3)
                        //
                        result[3 - 1] = zero;
                        for (i = 1; i <= n; i = i + 1) {
                            if (s[i - 1] == zero) {
                                if (dtru[i - 1] > abnorm * ulp) {
                                    result[3 - 1] = ulpinv;
                                }
                            } else if (dtru[i - 1] == zero) {
                                if (s[i - 1] > abnorm * ulp) {
                                    result[3 - 1] = ulpinv;
                                }
                            } else {
                                rwork[i - 1] = max(abs(dtru[i - 1] / s[i - 1]), abs(s[i - 1] / dtru[i - 1]));
                                result[3 - 1] = max(result[3 - 1], rwork[i - 1]);
                            }
                        }
                        //
                        //                    Test (4)
                        //
                        result[4 - 1] = zero;
                        if (dif[1 - 1] == zero) {
                            if (diftru[1 - 1] > abnorm * ulp) {
                                result[4 - 1] = ulpinv;
                            }
                        } else if (diftru[1 - 1] == zero) {
                            if (dif[1 - 1] > abnorm * ulp) {
                                result[4 - 1] = ulpinv;
                            }
                        } else if (dif[5 - 1] == zero) {
                            if (diftru[5 - 1] > abnorm * ulp) {
                                result[4 - 1] = ulpinv;
                            }
                        } else if (diftru[5 - 1] == zero) {
                            if (dif[5 - 1] > abnorm * ulp) {
                                result[4 - 1] = ulpinv;
                            }
                        } else {
                            ratio1 = max(abs(diftru[1 - 1] / dif[1 - 1]), abs(dif[1 - 1] / diftru[1 - 1]));
                            ratio2 = max(abs(diftru[5 - 1] / dif[5 - 1]), abs(dif[5 - 1] / diftru[5 - 1]));
                            result[4 - 1] = max(ratio1, ratio2);
                        }
                        //
                        ntestt += 4;
                        //
                        //                    Print out tests which fail.
                        //
                        for (j = 1; j <= 4; j = j + 1) {
                            if ((result[j - 1] >= thrsh2 && j >= 4) || (result[j - 1] >= thresh && j <= 3)) {
                                //
                                //                       If this is the first test to fail,
                                //                       prINTEGER a header to the data file.
                                //
                                if (nerrs == 0) {
                                    write(nout, format_9997), "ZXV";
                                    //
                                    //                          Print out messages for built-in examples
                                    //
                                    //                          Matrix types
                                    //
                                    write(nout, "(' Matrix types: ',/)");
                                    write(nout, "(' TYPE 1: Da is diagonal, Db is identity, ',/,"
                                                "'     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ',/,"
                                                "'     YH and X are left and right eigenvectors. ',/)");
                                    write(nout, "(' TYPE 2: Da is quasi-diagonal, Db is identity, ',/,"
                                                "'     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) ',/,"
                                                "'     YH and X are left and right eigenvectors. ',/)");
                                    //
                                    //                          Tests performed
                                    //
                                    write(nout, format_9992), "'", "transpose", "'";
                                    //
                                }
                                nerrs++;
                                if (result[j - 1] < 10000.0) {
                                    sprintnum_short(buf, result[j - 1]);
                                    write(nout, "(' Type=',i2,',',' IWA=',i2,', IWB=',i2,', IWX=',i2,"
                                                "', IWY=',i2,', result ',i2,' is',0p,buf)"),
                                        iptype, iwa, iwb, iwx, iwy, j, buf;
                                } else {
                                    sprintnum_short(buf, result[j - 1]);
                                    write(nout, "(' Type=',i2,',',' IWA=',i2,', IWB=',i2,', IWX=',i2,"
                                                "', IWY=',i2,', result ',i2,' is',1p,buf)"),
                                        iptype, iwa, iwb, iwx, iwy, j, buf;
                                }
                            }
                        }
                    //
                    statement_30:;
                        //
                    }
                }
            }
        }
    }
    //
    goto statement_150;
//
statement_90:
    //
    //     Read in data from file to check accuracy of condition estimation
    //     Read input data until N=0
    //
    try {
        read(nin, star), n;
    } catch (read_end const) {
        goto statement_150;
    }
    if (n == 0) {
        goto statement_150;
    }
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, a(i, j);
            }
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, b(i, j);
            }
        }
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= n; i = i + 1) {
            rloop, dtru(i);
        }
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= n; i = i + 1) {
            rloop, diftru(i);
        }
    }
    //
    nptknt++;
    //
    //     Compute eigenvalues/eigenvectors of (A, B).
    //     Compute eigenvalue/eigenvector condition numbers
    //     using computed eigenvectors.
    //
    Clacpy("F", n, n, a, lda, ai, lda);
    Clacpy("F", n, n, b, lda, bi, lda);
    //
    Cggevx("N", "V", "V", "B", n, ai, lda, bi, lda, alpha, beta, vl, lda, vr, lda, ilo, ihi, lscale, rscale, anorm, bnorm, s, dif, work, lwork, rwork, iwork, bwork, linfo);
    //
    if (linfo != 0) {
        write(nout, "(' Cdrgvx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,"
                    "', Input example #',i2,')')"),
            "Cggevx", linfo, n, nptknt;
        goto statement_140;
    }
    //
    //     Compute the norm(A, B)
    //
    Clacpy("Full", n, n, ai, lda, work, n);
    Clacpy("Full", n, n, bi, lda, &work[(n * n + 1) - 1], n);
    abnorm = Clange("Fro", n, 2 * n, work, n, rwork);
    //
    //     Tests (1) and (2)
    //
    result[1 - 1] = zero;
    Cget52(true, n, a, lda, b, lda, vl, lda, alpha, beta, work, rwork, result[1 - 1]);
    if (result[2 - 1] > thresh) {
        write(nout, format_9986), "Left", "Cggevx", result(2), n, nptknt;
    }
    //
    result[2 - 1] = zero;
    Cget52(false, n, a, lda, b, lda, vr, lda, alpha, beta, work, rwork, result[2 - 1]);
    if (result[3 - 1] > thresh) {
        write(nout, format_9986), "Right", "Cggevx", result(3), n, nptknt;
    }
    //
    //     Test (3)
    //
    result[3 - 1] = zero;
    for (i = 1; i <= n; i = i + 1) {
        if (s[i - 1] == zero) {
            if (dtru[i - 1] > abnorm * ulp) {
                result[3 - 1] = ulpinv;
            }
        } else if (dtru[i - 1] == zero) {
            if (s[i - 1] > abnorm * ulp) {
                result[3 - 1] = ulpinv;
            }
        } else {
            rwork[i - 1] = max(abs(dtru[i - 1] / s[i - 1]), abs(s[i - 1] / dtru[i - 1]));
            result[3 - 1] = max(result[3 - 1], &rwork[i - 1]);
        }
    }
    //
    //     Test (4)
    //
    result[4 - 1] = zero;
    if (dif[1 - 1] == zero) {
        if (diftru[1 - 1] > abnorm * ulp) {
            result[4 - 1] = ulpinv;
        }
    } else if (diftru[1 - 1] == zero) {
        if (dif[1 - 1] > abnorm * ulp) {
            result[4 - 1] = ulpinv;
        }
    } else if (dif[5 - 1] == zero) {
        if (diftru[5 - 1] > abnorm * ulp) {
            result[4 - 1] = ulpinv;
        }
    } else if (diftru[5 - 1] == zero) {
        if (dif[5 - 1] > abnorm * ulp) {
            result[4 - 1] = ulpinv;
        }
    } else {
        ratio1 = max(abs(diftru[1 - 1] / dif[1 - 1]), abs(dif[1 - 1] / diftru[1 - 1]));
        ratio2 = max(abs(diftru[5 - 1] / dif[5 - 1]), abs(dif[5 - 1] / diftru[5 - 1]));
        result[4 - 1] = max(ratio1, ratio2);
    }
    //
    ntestt += 4;
    //
    //     Print out tests which fail.
    //
    for (j = 1; j <= 4; j = j + 1) {
        if (result[j - 1] >= thrsh2) {
            //
            //           If this is the first test to fail,
            //           prINTEGER a header to the data file.
            //
            if (nerrs == 0) {
                write(nout, format_9997), "ZXV";
                //
                //              Print out messages for built-in examples
                //
                //              Matrix types
                //
                write(nout, "('Input Example')");
                //
                //              Tests performed
                //
                write(nout, format_9992), "'", "transpose", "'";
                //
            }
            nerrs++;
            if (result[j - 1] < 10000.0) {
                write(nout, "(' Input example #',i2,', matrix order=',i4,',',' result ',i2,"
                            "' is',0p,f8.2)"),
                    nptknt, n, j, result(j);
            } else {
                write(nout, "(' Input example #',i2,', matrix order=',i4,',',' result ',i2,"
                            "' is',1p,d10.3)"),
                    nptknt, n, j, result(j);
            }
        }
    }
//
statement_140:
    //
    goto statement_90;
statement_150:
    //
    //     Summary
    //
    Alasvm("ZXV", nout, nerrs, ntestt, 0);
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Cdrgvx
    //
}

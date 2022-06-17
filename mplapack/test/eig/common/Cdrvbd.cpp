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

void Cdrvbd(INTEGER const nsizes, INTEGER *mm, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, COMPLEX *a, INTEGER const lda, COMPLEX *u, INTEGER const ldu, COMPLEX *vt, INTEGER const ldvt, COMPLEX *asav, COMPLEX *usav, COMPLEX *vtsav, REAL *s, REAL *ssav, REAL *e, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER const nounit, INTEGER &info) {
    INTEGER ldasav = lda;
    INTEGER ldusav = ldu;
    INTEGER ldvtsav = ldvt;

    common cmn;
    common_write write(cmn);
    char cjob[4] = {'N', 'O', 'S', 'A'};
    char cjobr[3] = {'A', 'V', 'I'};
    char cjobv[2] = {'N', 'V'};
    char buf[1024];
    INTEGER nerrs = 0;
    INTEGER ntestt = 0;
    INTEGER ntestf = 0;
    bool badmm = false;
    bool badnn = false;
    INTEGER mmax = 0;
    INTEGER nmax = 0;
    INTEGER mnmax = 0;
    INTEGER minwrk = 0;
    INTEGER j = 0;
    REAL unfl = 0.0;
    const REAL one = 1.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    REAL ulpinv = 0.0;
    REAL rtunfl = 0.0;
    INTEGER jsize = 0;
    INTEGER m = 0;
    INTEGER n = 0;
    INTEGER mnmin = 0;
    const INTEGER maxtyp = 5;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ntest = 0;
    INTEGER ioldsd[4];
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER i = 0;
    const REAL zero = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    REAL anorm = 0.0;
    INTEGER iinfo = 0;
    INTEGER iwspc = 0;
    INTEGER iwtmp = 0;
    INTEGER lswork = 0;
    REAL result[39];
    INTEGER iju = 0;
    INTEGER ijvt = 0;
    char jobu;
    char jobvt;
    REAL dif = 0.0;
    REAL div = 0.0;
    INTEGER ijq = 0;
    char jobq;
    INTEGER lrwork = 0;
    INTEGER liwork = 0;
    INTEGER numrank = 0;
    REAL vl = 0.0;
    REAL vu = 0.0;
    INTEGER il = 0;
    INTEGER iu = 0;
    INTEGER ns = 0;
    char range;
    INTEGER iseed2[4];
    INTEGER itemp = 0;
    INTEGER nsi = 0;
    const REAL half = 0.5e0;
    const REAL two = 2.0;
    INTEGER nsv = 0;
    INTEGER nfail = 0;
    static const char *format_9995 = "(' Cdrvbd: ',a,' returned INFO=',i6,'.',/,9x,'M=',i6,', N=',i6,"
                                     "', JTYPE=',i6,', LSWORK=',i6,/,9x,'ISEED=(',3(i5,','),i5,')')";
    //
    //     Check for errors
    //
    info = 0;
    //
    //     Important constants
    //
    nerrs = 0;
    ntestt = 0;
    ntestf = 0;
    badmm = false;
    badnn = false;
    mmax = 1;
    nmax = 1;
    mnmax = 1;
    minwrk = 1;
    for (j = 1; j <= nsizes; j = j + 1) {
        mmax = max(mmax, mm[j - 1]);
        if (mm[j - 1] < 0) {
            badmm = true;
        }
        nmax = max(nmax, nn[j - 1]);
        if (nn[j - 1] < 0) {
            badnn = true;
        }
        mnmax = max({mnmax, min(mm[j - 1], nn[j - 1])});
        minwrk = max({minwrk, max({3 * min(mm[j - 1], nn[j - 1]) + max(mm[j - 1], nn[j - 1]) * max(mm[j - 1], nn[j - 1]), 5 * min(mm[j - 1], nn[j - 1]), 3 * max(mm[j - 1], nn[j - 1])})});
    }
    //
    //     Check for errors
    //
    if (nsizes < 0) {
        info = -1;
    } else if (badmm) {
        info = -2;
    } else if (badnn) {
        info = -3;
    } else if (ntypes < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, mmax)) {
        info = -10;
    } else if (ldu < max((INTEGER)1, mmax)) {
        info = -12;
    } else if (ldvt < max((INTEGER)1, nmax)) {
        info = -14;
    } else if (minwrk > lwork) {
        info = -21;
    }
    //
    if (info != 0) {
        Mxerbla("Cdrvbd", -info);
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
    unfl = Rlamch("S");
    ovfl = one / unfl;
    ulp = Rlamch("E");
    ulpinv = one / ulp;
    rtunfl = sqrt(unfl);
    //
    //     Loop over sizes, types
    //
    nerrs = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        m = mm[jsize - 1];
        n = nn[jsize - 1];
        mnmin = min(m, n);
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_220;
            }
            ntest = 0;
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            //           Compute "A"
            //
            if (mtypes > maxtyp) {
                goto statement_50;
            }
            //
            if (jtype == 1) {
                //
                //              Zero matrix
                //
                Claset("Full", m, n, czero, czero, a, lda);
                for (i = 1; i <= min(m, n); i = i + 1) {
                    s[i - 1] = zero;
                }
                //
            } else if (jtype == 2) {
                //
                //              Identity matrix
                //
                Claset("Full", m, n, czero, cone, a, lda);
                for (i = 1; i <= min(m, n); i = i + 1) {
                    s[i - 1] = one;
                }
                //
            } else {
                //
                //              (Scaled) random matrix
                //
                if (jtype == 3) {
                    anorm = one;
                }
                if (jtype == 4) {
                    anorm = unfl / ulp;
                }
                if (jtype == 5) {
                    anorm = ovfl * ulp;
                }
                Clatms(m, n, "U", iseed, "N", s, 4, castREAL(mnmin), anorm, m - 1, n - 1, "N", a, lda, work, iinfo);
                if (iinfo != 0) {
                    write(nounit, "(' Cdrvbd: ',a,' returned INFO=',i6,'.',/,9x,'M=',i6,', N=',i6,"
                                  "', JTYPE=',i6,', ISEED=(',3(i5,','),i5,')')"),
                        "Generator", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
            }
        //
        statement_50:
            Clacpy("F", m, n, a, lda, asav, lda);
            //
            //           Do for minimal and adequate (for blocking) workspace
            //
            for (iwspc = 1; iwspc <= 4; iwspc = iwspc + 1) {
                //
                //              Test for Cgesvd
                //
                iwtmp = 2 * min(m, n) + max(m, n);
                lswork = iwtmp + (iwspc - 1) * (lwork - iwtmp) / 3;
                lswork = min(lswork, lwork);
                lswork = max(lswork, (INTEGER)1);
                if (iwspc == 4) {
                    lswork = lwork;
                }
                //
                for (j = 1; j <= 35; j = j + 1) {
                    result[j - 1] = -one;
                }
                //
                //              Factorize A
                //
                if (iwspc > 1) {
                    Clacpy("F", m, n, asav, lda, a, lda);
                }
                Cgesvd("A", "A", m, n, a, lda, ssav, usav, ldu, vtsav, ldvt, work, lswork, rwork, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9995), "GESVD", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                //              Do tests 1--4
                //
                Cbdt01(m, n, 0, asav, lda, usav, ldu, ssav, e, vtsav, ldvt, work, rwork, result[1 - 1]);
                if (m != 0 && n != 0) {
                    Cunt01("Columns", mnmin, m, usav, ldu, work, lwork, rwork, result[2 - 1]);
                    Cunt01("Rows", mnmin, n, vtsav, ldvt, work, lwork, rwork, result[3 - 1]);
                }
                result[4 - 1] = 0.0;
                for (i = 1; i <= mnmin - 1; i = i + 1) {
                    if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                        result[4 - 1] = ulpinv;
                    }
                    if (ssav[i - 1] < zero) {
                        result[4 - 1] = ulpinv;
                    }
                }
                if (mnmin >= 1) {
                    if (ssav[mnmin - 1] < zero) {
                        result[4 - 1] = ulpinv;
                    }
                }
                //
                //              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
                //
                result[5 - 1] = zero;
                result[6 - 1] = zero;
                result[7 - 1] = zero;
                for (iju = 0; iju <= 3; iju = iju + 1) {
                    for (ijvt = 0; ijvt <= 3; ijvt = ijvt + 1) {
                        if ((iju == 3 && ijvt == 3) || (iju == 1 && ijvt == 1)) {
                            goto statement_90;
                        }
                        jobu = cjob[(iju + 1) - 1];
                        jobvt = cjob[(ijvt + 1) - 1];
                        Clacpy("F", m, n, asav, lda, a, lda);
                        Cgesvd(&jobu, &jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lswork, rwork, iinfo);
                        //
                        //                    Compare U
                        //
                        dif = zero;
                        if (m > 0 && n > 0) {
                            if (iju == 1) {
                                Cunt03("C", m, mnmin, m, mnmin, usav, ldu, a, lda, work, lwork, rwork, dif, iinfo);
                            } else if (iju == 2) {
                                Cunt03("C", m, mnmin, m, mnmin, usav, ldu, u, ldu, work, lwork, rwork, dif, iinfo);
                            } else if (iju == 3) {
                                Cunt03("C", m, m, m, mnmin, usav, ldu, u, ldu, work, lwork, rwork, dif, iinfo);
                            }
                        }
                        result[5 - 1] = max(result[5 - 1], dif);
                        //
                        //                    Compare VT
                        //
                        dif = zero;
                        if (m > 0 && n > 0) {
                            if (ijvt == 1) {
                                Cunt03("R", n, mnmin, n, mnmin, vtsav, ldvt, a, lda, work, lwork, rwork, dif, iinfo);
                            } else if (ijvt == 2) {
                                Cunt03("R", n, mnmin, n, mnmin, vtsav, ldvt, vt, ldvt, work, lwork, rwork, dif, iinfo);
                            } else if (ijvt == 3) {
                                Cunt03("R", n, n, n, mnmin, vtsav, ldvt, vt, ldvt, work, lwork, rwork, dif, iinfo);
                            }
                        }
                        result[6 - 1] = max(result[6 - 1], dif);
                        //
                        //                    Compare S
                        //
                        dif = zero;
                        div = max(REAL(castREAL(mnmin) * ulp * s[1 - 1]), Rlamch("Safe minimum"));
                        for (i = 1; i <= mnmin - 1; i = i + 1) {
                            if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                                dif = ulpinv;
                            }
                            if (ssav[i - 1] < zero) {
                                dif = ulpinv;
                            }
                            dif = max(dif, REAL(abs(ssav[i - 1] - s[i - 1]) / div));
                        }
                        result[7 - 1] = max(result[7 - 1], dif);
                    statement_90:;
                    }
                }
                //
                //              Test for Cgesdd
                //
                iwtmp = 2 * mnmin * mnmin + 2 * mnmin + max(m, n);
                lswork = iwtmp + (iwspc - 1) * (lwork - iwtmp) / 3;
                lswork = min(lswork, lwork);
                lswork = max(lswork, (INTEGER)1);
                if (iwspc == 4) {
                    lswork = lwork;
                }
                //
                //              Factorize A
                //
                Clacpy("F", m, n, asav, lda, a, lda);
                Cgesdd("A", m, n, a, lda, ssav, usav, ldu, vtsav, ldvt, work, lswork, rwork, iwork, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9995), "GESDD", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                //              Do tests 1--4
                //
                Cbdt01(m, n, 0, asav, lda, usav, ldu, ssav, e, vtsav, ldvt, work, rwork, result[8 - 1]);
                if (m != 0 && n != 0) {
                    Cunt01("Columns", mnmin, m, usav, ldu, work, lwork, rwork, result[9 - 1]);
                    Cunt01("Rows", mnmin, n, vtsav, ldvt, work, lwork, rwork, result[10 - 1]);
                }
                result[11 - 1] = 0.0;
                for (i = 1; i <= mnmin - 1; i = i + 1) {
                    if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                        result[11 - 1] = ulpinv;
                    }
                    if (ssav[i - 1] < zero) {
                        result[11 - 1] = ulpinv;
                    }
                }
                if (mnmin >= 1) {
                    if (ssav[mnmin - 1] < zero) {
                        result[11 - 1] = ulpinv;
                    }
                }
                //
                //              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
                //
                result[12 - 1] = zero;
                result[13 - 1] = zero;
                result[14 - 1] = zero;
                for (ijq = 0; ijq <= 2; ijq = ijq + 1) {
                    jobq = cjob[(ijq + 1) - 1];
                    Clacpy("F", m, n, asav, lda, a, lda);
                    Cgesdd(&jobq, m, n, a, lda, s, u, ldu, vt, ldvt, work, lswork, rwork, iwork, iinfo);
                    //
                    //                 Compare U
                    //
                    dif = zero;
                    if (m > 0 && n > 0) {
                        if (ijq == 1) {
                            if (m >= n) {
                                Cunt03("C", m, mnmin, m, mnmin, usav, ldu, a, lda, work, lwork, rwork, dif, iinfo);
                            } else {
                                Cunt03("C", m, mnmin, m, mnmin, usav, ldu, u, ldu, work, lwork, rwork, dif, iinfo);
                            }
                        } else if (ijq == 2) {
                            Cunt03("C", m, mnmin, m, mnmin, usav, ldu, u, ldu, work, lwork, rwork, dif, iinfo);
                        }
                    }
                    result[12 - 1] = max(result[12 - 1], dif);
                    //
                    //                 Compare VT
                    //
                    dif = zero;
                    if (m > 0 && n > 0) {
                        if (ijq == 1) {
                            if (m >= n) {
                                Cunt03("R", n, mnmin, n, mnmin, vtsav, ldvt, vt, ldvt, work, lwork, rwork, dif, iinfo);
                            } else {
                                Cunt03("R", n, mnmin, n, mnmin, vtsav, ldvt, a, lda, work, lwork, rwork, dif, iinfo);
                            }
                        } else if (ijq == 2) {
                            Cunt03("R", n, mnmin, n, mnmin, vtsav, ldvt, vt, ldvt, work, lwork, rwork, dif, iinfo);
                        }
                    }
                    result[13 - 1] = max(result[13 - 1], dif);
                    //
                    //                 Compare S
                    //
                    dif = zero;
                    div = max(REAL(castREAL(mnmin) * ulp * s[1 - 1]), Rlamch("Safe minimum"));
                    for (i = 1; i <= mnmin - 1; i = i + 1) {
                        if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                            dif = ulpinv;
                        }
                        if (ssav[i - 1] < zero) {
                            dif = ulpinv;
                        }
                        dif = max(dif, REAL(abs(ssav[i - 1] - s[i - 1]) / div));
                    }
                    result[14 - 1] = max(result[14 - 1], dif);
                }
                //
                //              Test Cgesvdq
                //              Note: Cgesvdq only works for M >= N
                //
                result[36 - 1] = zero;
                result[37 - 1] = zero;
                result[38 - 1] = zero;
                result[39 - 1] = zero;
                //
                if (m >= n) {
                    iwtmp = 2 * mnmin * mnmin + 2 * mnmin + max(m, n);
                    lswork = iwtmp + (iwspc - 1) * (lwork - iwtmp) / 3;
                    lswork = min(lswork, lwork);
                    lswork = max(lswork, (INTEGER)1);
                    if (iwspc == 4) {
                        lswork = lwork;
                    }
                    //
                    Clacpy("F", m, n, asav, lda, a, lda);
                    //
                    lrwork = max({(INTEGER)2, m, 5 * n});
                    liwork = max(n, (INTEGER)1);
                    Cgesvdq("H", "N", "N", "A", "A", m, n, a, lda, ssav, usav, ldu, vtsav, ldvt, numrank, iwork, liwork, work, lwork, rwork, lrwork, iinfo);
                    //
                    if (iinfo != 0) {
                        write(nounit, format_9995), "Cgesvdq", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                        info = abs(iinfo);
                        return;
                    }
                    //
                    //                 Do tests 36--39
                    //
                    Cbdt01(m, n, 0, asav, lda, usav, ldu, ssav, e, vtsav, ldvt, work, rwork, result[36 - 1]);
                    if (m != 0 && n != 0) {
                        Cunt01("Columns", m, m, usav, ldu, work, lwork, rwork, result[37 - 1]);
                        Cunt01("Rows", n, n, vtsav, ldvt, work, lwork, rwork, result[38 - 1]);
                    }
                    result[39 - 1] = zero;
                    for (i = 1; i <= mnmin - 1; i = i + 1) {
                        if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                            result[39 - 1] = ulpinv;
                        }
                        if (ssav[i - 1] < zero) {
                            result[39 - 1] = ulpinv;
                        }
                    }
                    if (mnmin >= 1) {
                        if (ssav[mnmin - 1] < zero) {
                            result[39 - 1] = ulpinv;
                        }
                    }
                }
                //
                //              Test Cgesvj
                //              Note: Cgesvj only works for M >= N
                //
                result[15 - 1] = zero;
                result[16 - 1] = zero;
                result[17 - 1] = zero;
                result[18 - 1] = zero;
                //
                if (m >= n) {
                    iwtmp = 2 * mnmin * mnmin + 2 * mnmin + max(m, n);
                    lswork = iwtmp + (iwspc - 1) * (lwork - iwtmp) / 3;
                    lswork = min(lswork, lwork);
                    lswork = max(lswork, (INTEGER)1);
                    lrwork = max((INTEGER)6, n);
                    if (iwspc == 4) {
                        lswork = lwork;
                    }
                    //
                    Clacpy("F", m, n, asav, lda, usav, lda);
                    Cgesvj("G", "U", "V", m, n, usav, lda, ssav, 0, a, ldvt, work, lwork, rwork, lrwork, iinfo);
                    //
                    //                 Cgesvj returns V not VH
                    //
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= n; i = i + 1) {
                            vtsav[(j - 1) + (i - 1) * ldvtsav] = conj(a[(i - 1) + (j - 1) * lda]);
                        }
                    }
                    //
                    if (iinfo != 0) {
                        write(nounit, format_9995), "GESVJ", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                        info = abs(iinfo);
                        return;
                    }
                    //
                    //                 Do tests 15--18
                    //
                    Cbdt01(m, n, 0, asav, lda, usav, ldu, ssav, e, vtsav, ldvt, work, rwork, result[15 - 1]);
                    if (m != 0 && n != 0) {
                        Cunt01("Columns", m, m, usav, ldu, work, lwork, rwork, result[16 - 1]);
                        Cunt01("Rows", n, n, vtsav, ldvt, work, lwork, rwork, result[17 - 1]);
                    }
                    result[18 - 1] = zero;
                    for (i = 1; i <= mnmin - 1; i = i + 1) {
                        if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                            result[18 - 1] = ulpinv;
                        }
                        if (ssav[i - 1] < zero) {
                            result[18 - 1] = ulpinv;
                        }
                    }
                    if (mnmin >= 1) {
                        if (ssav[mnmin - 1] < zero) {
                            result[18 - 1] = ulpinv;
                        }
                    }
                }
                //
                //              Test Cgejsv
                //              Note: Cgejsv only works for M >= N
                //
                result[19 - 1] = zero;
                result[20 - 1] = zero;
                result[21 - 1] = zero;
                result[22 - 1] = zero;
                if (m >= n) {
                    iwtmp = 2 * mnmin * mnmin + 2 * mnmin + max(m, n);
                    lswork = iwtmp + (iwspc - 1) * (lwork - iwtmp) / 3;
                    lswork = min(lswork, lwork);
                    lswork = max(lswork, (INTEGER)1);
                    if (iwspc == 4) {
                        lswork = lwork;
                    }
                    lrwork = max((INTEGER)7, n + 2 * m);
                    //
                    Clacpy("F", m, n, asav, lda, vtsav, lda);
                    Cgejsv("G", "U", "V", "R", "N", "N", m, n, vtsav, lda, ssav, usav, ldu, a, ldvt, work, lwork, rwork, lrwork, iwork, iinfo);
                    //
                    //                 Cgejsv returns V not VH
                    //
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= n; i = i + 1) {
                            vtsav[(j - 1) + (i - 1) * ldvtsav] = conj(a[(i - 1) + (j - 1) * lda]);
                        }
                    }
                    //
                    if (iinfo != 0) {
                        write(nounit, format_9995), "GEJSV", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                        info = abs(iinfo);
                        return;
                    }
                    //
                    //                 Do tests 19--22
                    //
                    Cbdt01(m, n, 0, asav, lda, usav, ldu, ssav, e, vtsav, ldvt, work, rwork, result[19 - 1]);
                    if (m != 0 && n != 0) {
                        Cunt01("Columns", m, m, usav, ldu, work, lwork, rwork, result[20 - 1]);
                        Cunt01("Rows", n, n, vtsav, ldvt, work, lwork, rwork, result[21 - 1]);
                    }
                    result[22 - 1] = zero;
                    for (i = 1; i <= mnmin - 1; i = i + 1) {
                        if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                            result[22 - 1] = ulpinv;
                        }
                        if (ssav[i - 1] < zero) {
                            result[22 - 1] = ulpinv;
                        }
                    }
                    if (mnmin >= 1) {
                        if (ssav[mnmin - 1] < zero) {
                            result[22 - 1] = ulpinv;
                        }
                    }
                }
                //
                //              Test Cgesvdx
                //
                //              Factorize A
                //
                Clacpy("F", m, n, asav, lda, a, lda);
                Cgesvdx("V", "V", "A", m, n, a, lda, vl, vu, il, iu, ns, ssav, usav, ldu, vtsav, ldvt, work, lwork, rwork, iwork, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9995), "GESVDX", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                //              Do tests 1--4
                //
                result[23 - 1] = zero;
                result[24 - 1] = zero;
                result[25 - 1] = zero;
                Cbdt01(m, n, 0, asav, lda, usav, ldu, ssav, e, vtsav, ldvt, work, rwork, result[23 - 1]);
                if (m != 0 && n != 0) {
                    Cunt01("Columns", mnmin, m, usav, ldu, work, lwork, rwork, result[24 - 1]);
                    Cunt01("Rows", mnmin, n, vtsav, ldvt, work, lwork, rwork, result[25 - 1]);
                }
                result[26 - 1] = zero;
                for (i = 1; i <= mnmin - 1; i = i + 1) {
                    if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                        result[26 - 1] = ulpinv;
                    }
                    if (ssav[i - 1] < zero) {
                        result[26 - 1] = ulpinv;
                    }
                }
                if (mnmin >= 1) {
                    if (ssav[mnmin - 1] < zero) {
                        result[26 - 1] = ulpinv;
                    }
                }
                //
                //              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
                //
                result[27 - 1] = zero;
                result[28 - 1] = zero;
                result[29 - 1] = zero;
                for (iju = 0; iju <= 1; iju = iju + 1) {
                    for (ijvt = 0; ijvt <= 1; ijvt = ijvt + 1) {
                        if ((iju == 0 && ijvt == 0) || (iju == 1 && ijvt == 1)) {
                            goto statement_160;
                        }
                        jobu = cjobv[(iju + 1) - 1];
                        jobvt = cjobv[(ijvt + 1) - 1];
                        range = cjobr[1 - 1];
                        Clacpy("F", m, n, asav, lda, a, lda);
                        Cgesvdx(&jobu, &jobvt, "A", m, n, a, lda, vl, vu, il, iu, ns, ssav, u, ldu, vt, ldvt, work, lwork, rwork, iwork, iinfo);
                        //
                        //                    Compare U
                        //
                        dif = zero;
                        if (m > 0 && n > 0) {
                            if (iju == 1) {
                                Cunt03("C", m, mnmin, m, mnmin, usav, ldu, u, ldu, work, lwork, rwork, dif, iinfo);
                            }
                        }
                        result[27 - 1] = max(result[27 - 1], dif);
                        //
                        //                    Compare VT
                        //
                        dif = zero;
                        if (m > 0 && n > 0) {
                            if (ijvt == 1) {
                                Cunt03("R", n, mnmin, n, mnmin, vtsav, ldvt, vt, ldvt, work, lwork, rwork, dif, iinfo);
                            }
                        }
                        result[28 - 1] = max(result[28 - 1], dif);
                        //
                        //                    Compare S
                        //
                        dif = zero;
                        div = max(REAL(castREAL(mnmin) * ulp * s[1 - 1]), Rlamch("Safe minimum"));
                        for (i = 1; i <= mnmin - 1; i = i + 1) {
                            if (ssav[i - 1] < ssav[(i + 1) - 1]) {
                                dif = ulpinv;
                            }
                            if (ssav[i - 1] < zero) {
                                dif = ulpinv;
                            }
                            dif = max(dif, REAL(abs(ssav[i - 1] - s[i - 1]) / div));
                        }
                        result[29 - 1] = max(result[29 - 1], dif);
                    statement_160:;
                    }
                }
                //
                //              Do tests 8--10
                //
                for (i = 1; i <= 4; i = i + 1) {
                    iseed2[i - 1] = iseed[i - 1];
                }
                if (mnmin <= 1) {
                    il = 1;
                    iu = max((INTEGER)1, mnmin);
                } else {
                    il = 1 + castINTEGER((mnmin - 1) * Rlarnd(1, iseed2));
                    iu = 1 + castINTEGER((mnmin - 1) * Rlarnd(1, iseed2));
                    if (iu < il) {
                        itemp = iu;
                        iu = il;
                        il = itemp;
                    }
                }
                Clacpy("F", m, n, asav, lda, a, lda);
                Cgesvdx("V", "V", "I", m, n, a, lda, vl, vu, il, iu, nsi, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9995), "GESVDX", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                result[30 - 1] = zero;
                result[31 - 1] = zero;
                result[32 - 1] = zero;
                Cbdt05(m, n, asav, lda, s, nsi, u, ldu, vt, ldvt, work, result[30 - 1]);
                if (m != 0 && n != 0) {
                    Cunt01("Columns", m, nsi, u, ldu, work, lwork, rwork, result[31 - 1]);
                    Cunt01("Rows", nsi, n, vt, ldvt, work, lwork, rwork, result[32 - 1]);
                }
                //
                //              Do tests 11--13
                //
                if (mnmin > 0 && nsi > 1) {
                    if (il != 1) {
                        vu = ssav[il - 1] + max({REAL(half * abs(ssav[il - 1] - ssav[(il - 1) - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                    } else {
                        vu = ssav[1 - 1] + max({REAL(half * abs(ssav[ns - 1] - ssav[1 - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                    }
                    if (iu != ns) {
                        vl = ssav[iu - 1] - max({REAL(ulp * anorm), REAL(two * rtunfl), REAL(half * abs(ssav[(iu + 1) - 1] - ssav[iu - 1]))});
                    } else {
                        vl = ssav[ns - 1] - max({REAL(ulp * anorm), REAL(two * rtunfl), REAL(half * abs(ssav[ns - 1] - ssav[1 - 1]))});
                    }
                    vl = max(vl, zero);
                    vu = max(vu, zero);
                    if (vl >= vu) {
                        vu = max(REAL(vu * 2), REAL(vu + vl + half));
                    }
                } else {
                    vl = zero;
                    vu = one;
                }
                Clacpy("F", m, n, asav, lda, a, lda);
                Cgesvdx("V", "V", "V", m, n, a, lda, vl, vu, il, iu, nsv, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, iinfo);
                if (iinfo != 0) {
                    write(nounit, format_9995), "GESVDX", iinfo, m, n, jtype, lswork, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                result[33 - 1] = zero;
                result[34 - 1] = zero;
                result[35 - 1] = zero;
                Cbdt05(m, n, asav, lda, s, nsv, u, ldu, vt, ldvt, work, result[33 - 1]);
                if (m != 0 && n != 0) {
                    Cunt01("Columns", m, nsv, u, ldu, work, lwork, rwork, result[34 - 1]);
                    Cunt01("Rows", nsv, n, vt, ldvt, work, lwork, rwork, result[35 - 1]);
                }
                //
                //              End of Loop -- Check for RESULT(j) > THRESH
                //
                ntest = 0;
                nfail = 0;
                for (j = 1; j <= 39; j = j + 1) {
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
                    write(nounit, "(' SVD -- Complex Singular Value Decomposition Driver ',/,"
                                  "' Matrix types (see Cdrvbd for details):',/,/,' 1 = Zero matrix',"
                                  "/,' 2 = Identity matrix',/,"
                                  "' 3 = Evenly spaced singular values near 1',/,"
                                  "' 4 = Evenly spaced singular values near underflow',/,"
                                  "' 5 = Evenly spaced singular values near overflow',/,/,"
                                  "' Tests performed: ( A is dense, U and V are unitary,',/,19x,"
                                  "' S is an array, and Upartial, VTpartial, and',/,19x,"
                                  "' Spartial are partially computed U, VT and S),',/)");
                    sprintnum_short(buf, thresh);
                    write(nounit, "(' Tests performed with Test Threshold = ',a,/,' Cgesvd: ',/,"
                                  "' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/,"
                                  "' 2 = | I - U**T U | / ( M ulp ) ',/,"
                                  "' 3 = | I - VT VT**T | / ( N ulp ) ',/,"
                                  "' 4 = 0 if S contains min(M,N) nonnegative values in',"
                                  "' decreasing order, else 1/ulp',/,"
                                  "' 5 = | U - Upartial | / ( M ulp )',/,"
                                  "' 6 = | VT - VTpartial | / ( N ulp )',/,"
                                  "' 7 = | S - Spartial | / ( min(M,N) ulp |S| )',/,' Cgesdd: ',/,"
                                  "' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/,"
                                  "' 9 = | I - U**T U | / ( M ulp ) ',/,"
                                  "'10 = | I - VT VT**T | / ( N ulp ) ',/,"
                                  "'11 = 0 if S contains min(M,N) nonnegative values in',"
                                  "' decreasing order, else 1/ulp',/,"
                                  "'12 = | U - Upartial | / ( M ulp )',/,"
                                  "'13 = | VT - VTpartial | / ( N ulp )',/,"
                                  "'14 = | S - Spartial | / ( min(M,N) ulp |S| )',/,' Cgesvj: ',/,/,"
                                  "'15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/,"
                                  "'16 = | I - U**T U | / ( M ulp ) ',/,"
                                  "'17 = | I - VT VT**T | / ( N ulp ) ',/,"
                                  "'18 = 0 if S contains min(M,N) nonnegative values in',"
                                  "' decreasing order, else 1/ulp',/,' ZGESJV: ',/,/,"
                                  "'19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )',/,"
                                  "'20 = | I - U**T U | / ( M ulp ) ',/,"
                                  "'21 = | I - VT VT**T | / ( N ulp ) ',/,"
                                  "'22 = 0 if S contains min(M,N) nonnegative values in',"
                                  "' decreasing order, else 1/ulp',/,' Cgesvdx(V,V,A): ',/,"
                                  "'23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/,"
                                  "'24 = | I - U**T U | / ( M ulp ) ',/,"
                                  "'25 = | I - VT VT**T | / ( N ulp ) ',/,"
                                  "'26 = 0 if S contains min(M,N) nonnegative values in',"
                                  "' decreasing order, else 1/ulp',/,"
                                  "'27 = | U - Upartial | / ( M ulp )',/,"
                                  "'28 = | VT - VTpartial | / ( N ulp )',/,"
                                  "'29 = | S - Spartial | / ( min(M,N) ulp |S| )',/,"
                                  "' Cgesvdx(V,V,I): ',/,"
                                  "'30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/,"
                                  "'31 = | I - U**T U | / ( M ulp ) ',/,"
                                  "'32 = | I - VT VT**T | / ( N ulp ) ',/,' Cgesvdx(V,V,V) ',/,"
                                  "'33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/,"
                                  "'34 = | I - U**T U | / ( M ulp ) ',/,"
                                  "'35 = | I - VT VT**T | / ( N ulp ) ',' Cgesvdq(H,N,N,A,A',/,"
                                  "'36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/,"
                                  "'37 = | I - U**T U | / ( M ulp ) ',/,"
                                  "'38 = | I - VT VT**T | / ( N ulp ) ',/,"
                                  "'39 = 0 if S contains min(M,N) nonnegative values in',"
                                  "' decreasing order, else 1/ulp',/,/)"),
                        buf;
                    ntestf = 2;
                }
                //
                for (j = 1; j <= 39; j = j + 1) {
                    if (result[j - 1] >= thresh) {
                        sprintnum_short(buf, result[j - 1]);
                        write(nounit, "(' M=',i5,', N=',i5,', type ',i1,', IWS=',i1,', seed=',4(i4,"
                                      "','),' test(',i2,')=',a)"),
                            m, n, jtype, iwspc, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3], j, buf;
                    }
                }
                //
                nerrs += nfail;
                ntestt += ntest;
                //
            }
        //
        statement_220:;
        }
    }
    //
    //     Summary
    //
    Alasvm("ZBD", nounit, nerrs, ntestt, 0);
    //
    //     End of Cdrvbd
    //
}

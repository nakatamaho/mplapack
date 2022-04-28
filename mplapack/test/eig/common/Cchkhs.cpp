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

void Cchkhs(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, COMPLEX *a, INTEGER const lda, COMPLEX *h, COMPLEX *t1, COMPLEX *t2, COMPLEX *u, INTEGER const ldu, COMPLEX *z, COMPLEX *uz, COMPLEX *w1, COMPLEX *w3, COMPLEX *evectl, COMPLEX *evectr, COMPLEX *evecty, COMPLEX *evectx, COMPLEX *uu, COMPLEX *tau, COMPLEX *work, INTEGER const nwork, REAL *rwork, INTEGER *iwork, bool *select, REAL *result, INTEGER &info) {
    INTEGER ldh = lda;
    INTEGER ldt1 = lda;
    INTEGER ldt2 = lda;
    INTEGER ldz = ldu;
    INTEGER lduz = ldu;
    INTEGER ldevectl = ldu;
    INTEGER ldevectr = ldu;
    INTEGER ldevecty = ldu;
    INTEGER ldevectx = ldu;
    INTEGER lduu = ldu;
    char buf[1024];

    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 21;
    INTEGER ktype[21] = {1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9};
    INTEGER kmagn[21] = {1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3};
    INTEGER kmode[21] = {0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1};
    INTEGER kconds[21] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0};
    INTEGER ntestt = 0;
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL unfl = 0.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL ulpinv = 0.0;
    REAL rtunfl = 0.0;
    REAL rtovfl = 0.0;
    REAL rtulp = 0.0;
    REAL rtulpi = 0.0;
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER n1 = 0;
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
    REAL conds = 0.0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    INTEGER i = 0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    COMPLEX cdumma[4];
    INTEGER in = 0;
    REAL dumma[4];
    INTEGER k = 0;
    bool match = false;
    INTEGER jj = 0;
    static const char *format_9997 = "(' Cchkhs: Selected ',a,' Eigenvectors from ',a,"
                                     "' do not match other eigenvectors ',9x,'N=',i6,', JTYPE=',i6,', ISEED=(',"
                                     "3(i5,','),i5,')')";
    static const char *format_9998 = "(' Cchkhs: ',a,' Eigenvectors from ',a,' incorrectly ','normalized.',/,"
                                     "' Bits of error=',0p,a,',',9x,'N=',i6,', JTYPE=',i6,', ISEED=(',3(i5,"
                                     "','),i5,')')";
    static const char *format_9999 = "(' Cchkhs: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
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
    } else if (lda <= 1 || lda < nmax) {
        info = -9;
    } else if (ldu <= 1 || ldu < nmax) {
        info = -14;
    } else if (4 * nmax * nmax + 2 > nwork) {
        info = -26;
    }
    //
    if (info != 0) {
        Mxerbla("Cchkhs", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (nsizes == 0 || ntypes == 0) {
        return;
    }
    //
    //     More important constants
    //
    unfl = Rlamch("Safe minimum");
    ovfl = Rlamch("Overflow");
    Rlabad(unfl, ovfl);
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    ulpinv = one / ulp;
    rtunfl = sqrt(unfl);
    rtovfl = sqrt(ovfl);
    rtulp = sqrt(ulp);
    rtulpi = one / rtulp;
    //
    //     Loop over sizes, types
    //
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        if (n == 0) {
            goto statement_260;
        }
        n1 = max((INTEGER)1, n);
        aninv = one / castREAL(n1);
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_250;
            }
            nmats++;
            ntest = 0;
            //
            //           Save ISEED in case of an error.
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            //           Initialize RESULT
            //
            for (j = 1; j <= 14; j = j + 1) {
                result[j - 1] = zero;
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
            //       =5                 random log   hermitian, w/ eigenvalues
            //       =6                 random       general, w/ eigenvalues
            //       =7                              random diagonal
            //       =8                              random hermitian
            //       =9                              random general
            //       =10                             random triangular
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
            cond = ulpinv;
            //
            //           Special Matrices
            //
            if (itype == 1) {
                //
                //              Zero
                //
                iinfo = 0;
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
                Clatmr(n, n, "D", iseed, "N", work, imode, cond, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 5) {
                //
                //              Hermitian, eigenvalues specified
                //
                Clatms(n, n, "D", iseed, "H", rwork, imode, cond, anorm, n, n, "N", a, lda, work, iinfo);
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
                Clatme(n, "D", iseed, work, imode, cond, cone, "T", "T", "T", rwork, 4, conds, n, n, anorm, a, lda, &work[(n + 1) - 1], iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random eigenvalues
                //
                Clatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Hermitian, random eigenvalues
                //
                Clatmr(n, n, "D", iseed, "H", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              General, random eigenvalues
                //
                Clatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 10) {
                //
                //              Triangular, random eigenvalues
                //
                Clatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
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
            //           Call Cgehrd to compute H and U, do tests.
            //
            Clacpy(" ", n, n, a, lda, h, lda);
            ntest = 1;
            //
            ilo = 1;
            ihi = n;
            //
            Cgehrd(n, ilo, ihi, h, lda, work, &work[(n + 1) - 1], nwork - n, iinfo);
            //
            if (iinfo != 0) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Cgehrd", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                goto statement_240;
            }
            //
            for (j = 1; j <= n - 1; j = j + 1) {
                uu[((j + 1) - 1) + (j - 1) * lduu] = czero;
                for (i = j + 2; i <= n; i = i + 1) {
                    u[(i - 1) + (j - 1) * ldu] = h[(i - 1) + (j - 1) * ldh];
                    uu[(i - 1) + (j - 1) * lduu] = h[(i - 1) + (j - 1) * ldh];
                    h[(i - 1) + (j - 1) * ldh] = czero;
                }
            }
            Ccopy(n - 1, work, 1, tau, 1);
            Cunghr(n, ilo, ihi, u, ldu, work, &work[(n + 1) - 1], nwork - n, iinfo);
            ntest = 2;
            //
            Chst01(n, ilo, ihi, a, lda, h, lda, u, ldu, work, nwork, rwork, &result[1 - 1]);
            //
            //           Call Chseqr to compute T1, T2 and Z, do tests.
            //
            //           Eigenvalues only (W3)
            //
            Clacpy(" ", n, n, h, lda, t2, lda);
            ntest = 3;
            result[3 - 1] = ulpinv;
            //
            Chseqr("E", "N", n, ilo, ihi, t2, lda, w3, uz, ldu, work, nwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Chseqr(E)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                if (iinfo <= n + 2) {
                    info = abs(iinfo);
                    goto statement_240;
                }
            }
            //
            //           Eigenvalues (W1) and Full Schur Form (T2)
            //
            Clacpy(" ", n, n, h, lda, t2, lda);
            //
            Chseqr("S", "N", n, ilo, ihi, t2, lda, w1, uz, ldu, work, nwork, iinfo);
            if (iinfo != 0 && iinfo <= n + 2) {
                write(nounit, format_9999), "Chseqr(S)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                goto statement_240;
            }
            //
            //           Eigenvalues (W1), Schur Form (T1), and Schur Vectors (UZ)
            //
            Clacpy(" ", n, n, h, lda, t1, lda);
            Clacpy(" ", n, n, u, ldu, uz, ldu);
            //
            Chseqr("S", "V", n, ilo, ihi, t1, lda, w1, uz, ldu, work, nwork, iinfo);
            if (iinfo != 0 && iinfo <= n + 2) {
                write(nounit, format_9999), "Chseqr(V)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                goto statement_240;
            }
            //
            //           Compute Z = U' UZ
            //
            Cgemm("C", "N", n, n, n, cone, u, ldu, uz, ldu, czero, z, ldu);
            ntest = 8;
            //
            //           Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
            //                and 4: | I - Z Z' | / ( n ulp )
            //
            Chst01(n, ilo, ihi, h, lda, t1, lda, z, ldu, work, nwork, rwork, &result[3 - 1]);
            //
            //           Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
            //                and 6: | I - UZ (UZ)' | / ( n ulp )
            //
            Chst01(n, ilo, ihi, a, lda, t1, lda, uz, ldu, work, nwork, rwork, &result[5 - 1]);
            //
            //           Do Test 7: | T2 - T1 | / ( |T| n ulp )
            //
            Cget10(n, n, t2, lda, t1, lda, work, rwork, result[7 - 1]);
            //
            //           Do Test 8: | W3 - W1 | / ( max(|W1|,|W3|) ulp )
            //
            temp1 = zero;
            temp2 = zero;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = max({temp1, REAL(abs(w1[j - 1])), REAL(abs(w3[j - 1]))});
                temp2 = max(temp2, REAL(abs(w1[j - 1] - w3[j - 1])));
            }
            //
            result[8 - 1] = temp2 / max(unfl, REAL(ulp * max(temp1, temp2)));
            //
            //           Compute the Left and Right Eigenvectors of T
            //
            //           Compute the Right eigenvector Matrix:
            //
            ntest = 9;
            result[9 - 1] = ulpinv;
            //
            //           Select every other eigenvector
            //
            for (j = 1; j <= n; j = j + 1) {
                select[j - 1] = false;
            }
            for (j = 1; j <= n; j = j + 2) {
                select[j - 1] = true;
            }
            Ctrevc("Right", "All", select, n, t1, lda, cdumma, ldu, evectr, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctrevc(R,A)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                goto statement_240;
            }
            //
            //           Test 9:  | TR - RW | / ( |T| |R| ulp )
            //
            Cget22("N", "N", "N", n, t1, lda, evectr, ldu, w1, work, rwork, &dumma[1 - 1]);
            result[9 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thresh) {
                sprintnum_short(buf, dumma[2 - 1]);
                write(nounit, format_9998), "Right", "Ctrevc", buf, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
            }
            //
            //           Compute selected right eigenvectors and confirm that
            //           they agree with previous right eigenvectors
            //
            Ctrevc("Right", "Some", select, n, t1, lda, cdumma, ldu, evectl, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctrevc(R,S)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                goto statement_240;
            }
            //
            k = 1;
            match = true;
            for (j = 1; j <= n; j = j + 1) {
                if (select[j - 1]) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (evectr[(jj - 1) + (j - 1) * ldevectr] != evectl[(jj - 1) + (k - 1) * ldevectl]) {
                            match = false;
                            goto statement_180;
                        }
                    }
                    k++;
                }
            }
        statement_180:
            if (!match) {
                write(nounit, format_9997), "Right", "Ctrevc", n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
            }
            //
            //           Compute the Left eigenvector Matrix:
            //
            ntest = 10;
            result[10 - 1] = ulpinv;
            Ctrevc("Left", "All", select, n, t1, lda, evectl, ldu, cdumma, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctrevc(L,A)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                goto statement_240;
            }
            //
            //           Test 10:  | LT - WL | / ( |T| |L| ulp )
            //
            Cget22("C", "N", "C", n, t1, lda, evectl, ldu, w1, work, rwork, &dumma[3 - 1]);
            result[10 - 1] = dumma[3 - 1];
            if (dumma[4 - 1] > thresh) {
                sprintnum_short(buf, result[4 - 1]);
                write(nounit, format_9998), "Left", "Ctrevc", buf, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
            }
            //
            //           Compute selected left eigenvectors and confirm that
            //           they agree with previous left eigenvectors
            //
            Ctrevc("Left", "Some", select, n, t1, lda, evectr, ldu, cdumma, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctrevc(L,S)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                goto statement_240;
            }
            //
            k = 1;
            match = true;
            for (j = 1; j <= n; j = j + 1) {
                if (select[j - 1]) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (evectl[(jj - 1) + (j - 1) * ldevectl] != evectr[(jj - 1) + (k - 1) * ldevectr]) {
                            match = false;
                            goto statement_210;
                        }
                    }
                    k++;
                }
            }
        statement_210:
            if (!match) {
                write(nounit, format_9997), "Left", "Ctrevc", n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
            }
            //
            //           Call Chsein for Right eigenvectors of H, do test 11
            //
            ntest = 11;
            result[11 - 1] = ulpinv;
            for (j = 1; j <= n; j = j + 1) {
                select[j - 1] = true;
            }
            //
            Chsein("Right", "Qr", "Ninitv", select, n, h, lda, w3, cdumma, ldu, evectx, ldu, n1, in, work, rwork, iwork, iwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Chsein(R)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_240;
                }
            } else {
                //
                //              Test 11:  | HX - XW | / ( |H| |X| ulp )
                //
                //                        (from inverse iteration)
                //
                Cget22("N", "N", "N", n, h, lda, evectx, ldu, w3, work, rwork, &dumma[1 - 1]);
                if (dumma[1 - 1] < ulpinv) {
                    result[11 - 1] = dumma[1 - 1] * aninv;
                }
                if (dumma[2 - 1] > thresh) {
                    sprintnum_short(buf, dumma[2 - 1]);
                    write(nounit, format_9998), "Right", "Chsein", buf, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                }
            }
            //
            //           Call Chsein for Left eigenvectors of H, do test 12
            //
            ntest = 12;
            result[12 - 1] = ulpinv;
            for (j = 1; j <= n; j = j + 1) {
                select[j - 1] = true;
            }
            //
            Chsein("Left", "Qr", "Ninitv", select, n, h, lda, w3, evecty, ldu, cdumma, ldu, n1, in, work, rwork, iwork, iwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Chsein(L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_240;
                }
            } else {
                //
                //              Test 12:  | YH - WY | / ( |H| |Y| ulp )
                //
                //                        (from inverse iteration)
                //
                Cget22("C", "N", "C", n, h, lda, evecty, ldu, w3, work, rwork, &dumma[3 - 1]);
                if (dumma[3 - 1] < ulpinv) {
                    result[12 - 1] = dumma[3 - 1] * aninv;
                }
                if (dumma[4 - 1] > thresh) {
                    sprintnum_short(buf, dumma[4 - 1]);
                    write(nounit, format_9998), "Left", "Chsein", buf, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                }
            }
            //
            //           Call Cunmhr for Right eigenvectors of A, do test 13
            //
            ntest = 13;
            result[13 - 1] = ulpinv;
            //
            Cunmhr("Left", "No transpose", n, n, ilo, ihi, uu, ldu, tau, evectx, ldu, work, nwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cunmhr(L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_240;
                }
            } else {
                //
                //              Test 13:  | AX - XW | / ( |A| |X| ulp )
                //
                //                        (from inverse iteration)
                //
                Cget22("N", "N", "N", n, a, lda, evectx, ldu, w3, work, rwork, &dumma[1 - 1]);
                if (dumma[1 - 1] < ulpinv) {
                    result[13 - 1] = dumma[1 - 1] * aninv;
                }
            }
            //
            //           Call Cunmhr for Left eigenvectors of A, do test 14
            //
            ntest = 14;
            result[14 - 1] = ulpinv;
            //
            Cunmhr("Left", "No transpose", n, n, ilo, ihi, uu, ldu, tau, evecty, ldu, work, nwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cunmhr(L)", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_240;
                }
            } else {
                //
                //              Test 14:  | YA - WY | / ( |A| |Y| ulp )
                //
                //                        (from inverse iteration)
                //
                Cget22("C", "N", "C", n, a, lda, evecty, ldu, w3, work, rwork, &dumma[3 - 1]);
                if (dumma[3 - 1] < ulpinv) {
                    result[14 - 1] = dumma[3 - 1] * aninv;
                }
            }
        //
        //           End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_240:
            //
            ntestt += ntest;
            Rlafts("ZHS", n, n, jtype, ntest, result, ioldsd, thresh, nounit, nerrs);
        //
        statement_250:;
        }
    statement_260:;
    }
    //
    //     Summary
    //
    Rlasum("ZHS", nounit, nerrs, ntestt);
    //
    //     End of Cchkhs
    //
}

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

void Rchkhs(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *h, REAL *t1, REAL *t2, REAL *u, INTEGER const ldu, REAL *z, REAL *uz, REAL *wr1, REAL *wi1, REAL *wr2, REAL *wi2, REAL *wr3, REAL *wi3, REAL *evectl, REAL *evectr, REAL *evecty, REAL *evectx, REAL *uu, REAL *tau, REAL *work, INTEGER const nwork, INTEGER *iwork, bool *select, REAL *result, INTEGER &info) {
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
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    INTEGER jcol = 0;
    REAL conds = 0.0;
    char adumma[1];
    INTEGER idumma[1];
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    INTEGER i = 0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    INTEGER nselc = 0;
    INTEGER nselr = 0;
    REAL dumma[6];
    INTEGER in = 0;
    INTEGER k = 0;
    bool match = false;
    INTEGER jj = 0;
    static const char *format_9997 = "(' Rchkhs: Selected ',a,' Eigenvectors from ',a,"
                                     "' do not match other eigenvectors ',9x,'N=',i6,', JTYPE=',i6,', ISEED=(',"
                                     "3(i5,','),i5,')')";
    static const char *format_9998 = "(' Rchkhs: ',a,' Eigenvectors from ',a,' incorrectly ','normalized.',/,"
                                     "' Bits of error=',0p,a,',',9x,'N=',i6,', JTYPE=',i6,', ISEED=(',3(i5,"
                                     "','),i5,')')";
    static const char *format_9999 = "(' Rchkhs: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
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
        info = -28;
    }
    //
    if (info != 0) {
        Mxerbla("Rchkhs", -info);
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
            goto statement_270;
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
                goto statement_260;
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
            //       =5                 random log   symmetric, w/ eigenvalues
            //       =6                 random       general, w/ eigenvalues
            //       =7                              random diagonal
            //       =8                              random symmetric
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
            Rlaset("Full", lda, n, zero, zero, a, lda);
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
                write(nounit, format_9999), "Generator", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                return;
            }
        //
        statement_100:
            //
            //           Call Rgehrd to compute H and U, do tests.
            //
            Rlacpy(" ", n, n, a, lda, h, lda);
            //
            ntest = 1;
            //
            ilo = 1;
            ihi = n;
            //
            Rgehrd(n, ilo, ihi, h, lda, work, &work[(n + 1) - 1], nwork - n, iinfo);
            //
            if (iinfo != 0) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Rgehrd", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_250;
            }
            //
            for (j = 1; j <= n - 1; j = j + 1) {
                uu[((j + 1) - 1) + (j - 1) * lduu] = zero;
                for (i = j + 2; i <= n; i = i + 1) {
                    u[(i - 1) + (j - 1) * ldu] = h[(i - 1) + (j - 1) * ldh];
                    uu[(i - 1) + (j - 1) * lduu] = h[(i - 1) + (j - 1) * ldh];
                    h[(i - 1) + (j - 1) * ldh] = zero;
                }
            }
            Rcopy(n - 1, work, 1, tau, 1);
            Rorghr(n, ilo, ihi, u, ldu, work, &work[(n + 1) - 1], nwork - n, iinfo);
            ntest = 2;
            //
            Rhst01(n, ilo, ihi, a, lda, h, lda, u, ldu, work, nwork, &result[1 - 1]);
            //
            //           Call Rhseqr to compute T1, T2 and Z, do tests.
            //
            //           Eigenvalues only (WR3,WI3)
            //
            Rlacpy(" ", n, n, h, lda, t2, lda);
            ntest = 3;
            result[3 - 1] = ulpinv;
            //
            Rhseqr("E", "N", n, ilo, ihi, t2, lda, wr3, wi3, uz, ldu, work, nwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rhseqr(E)", iinfo, n, jtype, ioldsd;
                if (iinfo <= n + 2) {
                    info = abs(iinfo);
                    goto statement_250;
                }
            }
            //
            //           Eigenvalues (WR2,WI2) and Full Schur Form (T2)
            //
            Rlacpy(" ", n, n, h, lda, t2, lda);
            //
            Rhseqr("S", "N", n, ilo, ihi, t2, lda, wr2, wi2, uz, ldu, work, nwork, iinfo);
            if (iinfo != 0 && iinfo <= n + 2) {
                write(nounit, format_9999), "Rhseqr(S)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_250;
            }
            //
            //           Eigenvalues (WR1,WI1), Schur Form (T1), and Schur vectors
            //           (UZ)
            //
            Rlacpy(" ", n, n, h, lda, t1, lda);
            Rlacpy(" ", n, n, u, ldu, uz, ldu);
            //
            Rhseqr("S", "V", n, ilo, ihi, t1, lda, wr1, wi1, uz, ldu, work, nwork, iinfo);
            if (iinfo != 0 && iinfo <= n + 2) {
                write(nounit, format_9999), "Rhseqr(V)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_250;
            }
            //
            //           Compute Z = U' UZ
            //
            Rgemm("T", "N", n, n, n, one, u, ldu, uz, ldu, zero, z, ldu);
            ntest = 8;
            //
            //           Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
            //                and 4: | I - Z Z' | / ( n ulp )
            //
            Rhst01(n, ilo, ihi, h, lda, t1, lda, z, ldu, work, nwork, &result[3 - 1]);
            //
            //           Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
            //                and 6: | I - UZ (UZ)' | / ( n ulp )
            //
            Rhst01(n, ilo, ihi, a, lda, t1, lda, uz, ldu, work, nwork, &result[5 - 1]);
            //
            //           Do Test 7: | T2 - T1 | / ( |T| n ulp )
            //
            Rget10(n, n, t2, lda, t1, lda, work, result[7 - 1]);
            //
            //           Do Test 8: | W2 - W1 | / ( max(|W1|,|W2|) ulp )
            //
            temp1 = zero;
            temp2 = zero;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = max({temp1, abs(wr1[j - 1]) + abs(wi1[j - 1]), abs(wr2[j - 1]) + abs(wi2[j - 1])});
                temp2 = max(temp2, abs(wr1[j - 1] - wr2[j - 1]) + abs(wi1[j - 1] - wi2[j - 1]));
            }
            //
            result[8 - 1] = temp2 / max({unfl, ulp * max(temp1, temp2)});
            //
            //           Compute the Left and Right Eigenvectors of T
            //
            //           Compute the Right eigenvector Matrix:
            //
            ntest = 9;
            result[9 - 1] = ulpinv;
            //
            //           Select last max(N/4,1) real, max(N/4,1) complex eigenvectors
            //
            nselc = 0;
            nselr = 0;
            j = n;
        statement_140:
            if (wi1[j - 1] == zero) {
                if (nselr < max(n / 4, 1)) {
                    nselr++;
                    select[j - 1] = true;
                } else {
                    select[j - 1] = false;
                }
                j = j - 1;
            } else {
                if (nselc < max(n / 4, 1)) {
                    nselc++;
                    select[j - 1] = true;
                    select[(j - 1) - 1] = false;
                } else {
                    select[j - 1] = false;
                    select[(j - 1) - 1] = false;
                }
                j = j - 2;
            }
            if (j > 0) {
                goto statement_140;
            }
            //
            Rtrevc("Right", "All", select, n, t1, lda, dumma, ldu, evectr, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtrevc(R,A)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_250;
            }
            //
            //           Test 9:  | TR - RW | / ( |T| |R| ulp )
            //
            Rget22("N", "N", "N", n, t1, lda, evectr, ldu, wr1, wi1, work, &dumma[1 - 1]);
            result[9 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thresh) {
                sprintnum_short(buf, dumma[2 - 1]);
                write(nounit, format_9998), "Right", "Rtrevc", buf, n, jtype, ioldsd;
            }
            //
            //           Compute selected right eigenvectors and confirm that
            //           they agree with previous right eigenvectors
            //
            Rtrevc("Right", "Some", select, n, t1, lda, dumma, ldu, evectl, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtrevc(R,S)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_250;
            }
            //
            k = 1;
            match = true;
            for (j = 1; j <= n; j = j + 1) {
                if (select[j - 1] && wi1[j - 1] == zero) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (evectr[(jj - 1) + (j - 1) * ldevectr] != evectl[(jj - 1) + (k - 1) * ldevectl]) {
                            match = false;
                            goto statement_180;
                        }
                    }
                    k++;
                } else if (select[j - 1] && wi1[j - 1] != zero) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (evectr[(jj - 1) + (j - 1) * ldevectr] != evectl[(jj - 1) + (k - 1) * ldevectl] || evectr[(jj - 1) + ((j + 1) - 1) * ldevectr] != evectl[(jj - 1) + ((k + 1) - 1) * ldevectl]) {
                            match = false;
                            goto statement_180;
                        }
                    }
                    k += 2;
                }
            }
        statement_180:
            if (!match) {
                write(nounit, format_9997), "Right", "Rtrevc", n, jtype, ioldsd;
            }
            //
            //           Compute the Left eigenvector Matrix:
            //
            ntest = 10;
            result[10 - 1] = ulpinv;
            Rtrevc("Left", "All", select, n, t1, lda, evectl, ldu, dumma, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtrevc(L,A)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_250;
            }
            //
            //           Test 10:  | LT - WL | / ( |T| |L| ulp )
            //
            Rget22("Trans", "N", "Conj", n, t1, lda, evectl, ldu, wr1, wi1, work, &dumma[3 - 1]);
            result[10 - 1] = dumma[3 - 1];
            if (dumma[4 - 1] > thresh) {
                sprintnum_short(buf, dumma[4 - 1]);
                write(nounit, format_9998), "Left", "Rtrevc", buf, n, jtype, ioldsd;
            }
            //
            //           Compute selected left eigenvectors and confirm that
            //           they agree with previous left eigenvectors
            //
            Rtrevc("Left", "Some", select, n, t1, lda, evectr, ldu, dumma, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtrevc(L,S)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_250;
            }
            //
            k = 1;
            match = true;
            for (j = 1; j <= n; j = j + 1) {
                if (select[j - 1] && wi1[j - 1] == zero) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (evectl[(jj - 1) + (j - 1) * ldevectl] != evectr[(jj - 1) + (k - 1) * ldevectr]) {
                            match = false;
                            goto statement_220;
                        }
                    }
                    k++;
                } else if (select[j - 1] && wi1[j - 1] != zero) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (evectl[(jj - 1) + (j - 1) * ldevectl] != evectr[(jj - 1) + (k - 1) * ldevectr] || evectl[(jj - 1) + ((j + 1) - 1) * ldevectl] != evectr[(jj - 1) + ((k + 1) - 1) * ldevectr]) {
                            match = false;
                            goto statement_220;
                        }
                    }
                    k += 2;
                }
            }
        statement_220:
            if (!match) {
                write(nounit, format_9997), "Left", "Rtrevc", n, jtype, ioldsd;
            }
            //
            //           Call Rhsein for Right eigenvectors of H, do test 11
            //
            ntest = 11;
            result[11 - 1] = ulpinv;
            for (j = 1; j <= n; j = j + 1) {
                select[j - 1] = true;
            }
            //
            Rhsein("Right", "Qr", "Ninitv", select, n, h, lda, wr3, wi3, dumma, ldu, evectx, ldu, n1, in, work, iwork, iwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rhsein(R)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_250;
                }
            } else {
                //
                //              Test 11:  | HX - XW | / ( |H| |X| ulp )
                //
                //                        (from inverse iteration)
                //
                Rget22("N", "N", "N", n, h, lda, evectx, ldu, wr3, wi3, work, &dumma[1 - 1]);
                if (dumma[1 - 1] < ulpinv) {
                    result[11 - 1] = dumma[1 - 1] * aninv;
                }
                if (dumma[2 - 1] > thresh) {
                    sprintnum_short(buf, dumma[2 - 1]);
                    write(nounit, format_9998), "Right", "Rhsein", buf, n, jtype, ioldsd;
                }
            }
            //
            //           Call Rhsein for Left eigenvectors of H, do test 12
            //
            ntest = 12;
            result[12 - 1] = ulpinv;
            for (j = 1; j <= n; j = j + 1) {
                select[j - 1] = true;
            }
            //
            Rhsein("Left", "Qr", "Ninitv", select, n, h, lda, wr3, wi3, evecty, ldu, dumma, ldu, n1, in, work, iwork, iwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rhsein(L)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_250;
                }
            } else {
                //
                //              Test 12:  | YH - WY | / ( |H| |Y| ulp )
                //
                //                        (from inverse iteration)
                //
                Rget22("C", "N", "C", n, h, lda, evecty, ldu, wr3, wi3, work, &dumma[3 - 1]);
                if (dumma[3 - 1] < ulpinv) {
                    result[12 - 1] = dumma[3 - 1] * aninv;
                }
                if (dumma[4 - 1] > thresh) {
                    sprintnum_short(buf, dumma[4 - 1]);
                    write(nounit, format_9998), "Left", "Rhsein", buf, n, jtype, ioldsd;
                }
            }
            //
            //           Call Rormhr for Right eigenvectors of A, do test 13
            //
            ntest = 13;
            result[13 - 1] = ulpinv;
            //
            Rormhr("Left", "No transpose", n, n, ilo, ihi, uu, ldu, tau, evectx, ldu, work, nwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rormhr(R)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_250;
                }
            } else {
                //
                //              Test 13:  | AX - XW | / ( |A| |X| ulp )
                //
                //                        (from inverse iteration)
                //
                Rget22("N", "N", "N", n, a, lda, evectx, ldu, wr3, wi3, work, &dumma[1 - 1]);
                if (dumma[1 - 1] < ulpinv) {
                    result[13 - 1] = dumma[1 - 1] * aninv;
                }
            }
            //
            //           Call Rormhr for Left eigenvectors of A, do test 14
            //
            ntest = 14;
            result[14 - 1] = ulpinv;
            //
            Rormhr("Left", "No transpose", n, n, ilo, ihi, uu, ldu, tau, evecty, ldu, work, nwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rormhr(L)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                if (iinfo < 0) {
                    goto statement_250;
                }
            } else {
                //
                //              Test 14:  | YA - WY | / ( |A| |Y| ulp )
                //
                //                        (from inverse iteration)
                //
                Rget22("C", "N", "C", n, a, lda, evecty, ldu, wr3, wi3, work, &dumma[3 - 1]);
                if (dumma[3 - 1] < ulpinv) {
                    result[14 - 1] = dumma[3 - 1] * aninv;
                }
            }
        //
        //           End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_250:
            //
            ntestt += ntest;
            Rlafts("DHS", n, n, jtype, ntest, result, ioldsd, thresh, nounit, nerrs);
        //
        statement_260:;
        }
    statement_270:;
    }
    //
    //     Summary
    //
    Rlasum("DHS", nounit, nerrs, ntestt);
    //
    //     End of Rchkhs
    //
}

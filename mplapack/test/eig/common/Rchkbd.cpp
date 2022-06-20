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

void Rchkbd(INTEGER const nsizes, INTEGER *mval, INTEGER *nval, INTEGER const ntypes, bool *dotype, INTEGER const nrhs, INTEGER *iseed, REAL const thresh, REAL *a, INTEGER const lda, REAL *bd, REAL *be, REAL *s1, REAL *s2, REAL *x, INTEGER const ldx, REAL *y, REAL *z, REAL *q, INTEGER const ldq, REAL *pt, INTEGER const ldpt, REAL *u, REAL *vt, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const nout, INTEGER &info) {

    INTEGER ldy = ldx;
    INTEGER ldz = ldx;
    INTEGER ldu = ldpt;
    INTEGER ldvt = ldpt;
    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 16;
    INTEGER ktype[16] = {1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9, 10};
    INTEGER kmagn[16] = {1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 0};
    INTEGER kmode[16] = {0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0};
    char buf[1024];
    bool badmm = false;
    bool badnn = false;
    INTEGER mmax = 0;
    INTEGER nmax = 0;
    INTEGER mnmax = 0;
    INTEGER minwrk = 0;
    INTEGER j = 0;
    char path[4];
    INTEGER nfail = 0;
    INTEGER ntest = 0;
    REAL unfl = 0.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL ulpinv = 0.0;
    const REAL two = 2.0;
    INTEGER log2ui = 0;
    REAL rtunfl = 0.0;
    REAL rtovfl = 0.0;
    REAL abstol = 0.0;
    INTEGER jsize = 0;
    INTEGER m = 0;
    INTEGER n = 0;
    INTEGER mnmin = 0;
    REAL amninv = 0.0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ioldsd[4];
    REAL result[40];
    char uplo;
    INTEGER itype = 0;
    INTEGER imode = 0;
    REAL anorm = 0.0;
    const REAL zero = 0.0;
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    bool bidiag = false;
    INTEGER jcol = 0;
    REAL temp1 = 0.0;
    INTEGER mq = 0;
    INTEGER i = 0;
    REAL temp2 = 0.0;
    const REAL half = 0.5e0;
    REAL dumma[1];
    REAL dum[1];
    INTEGER idum[1];
    INTEGER iwbs = 0;
    INTEGER iwbd = 0;
    INTEGER iwbe = 0;
    INTEGER iwbz = 0;
    INTEGER iwwork = 0;
    INTEGER mnmin2 = 0;
    INTEGER ns1 = 0;
    INTEGER ns2 = 0;
    INTEGER iseed2[4];
    INTEGER il = 0;
    INTEGER iu = 0;
    INTEGER itemp = 0;
    REAL vu = 0.0;
    REAL vl = 0.0;
    static const char *format_9998 = "(' Rchkbd: ',a,' returned INFO=',i6,'.',/,9x,'M=',i6,', N=',i6,"
                                     "', JTYPE=',i6,', ISEED=(',3(i5,','),i5,')')";
    //
    //     Check for errors
    //
    info = 0;
    //
    badmm = false;
    badnn = false;
    mmax = 1;
    nmax = 1;
    mnmax = 1;
    minwrk = 1;
    for (j = 1; j <= nsizes; j = j + 1) {
        mmax = max(mmax, mval[j - 1]);
        if (mval[j - 1] < 0) {
            badmm = true;
        }
        nmax = max(nmax, nval[j - 1]);
        if (nval[j - 1] < 0) {
            badnn = true;
        }
        mnmax = max({mnmax, min(mval[j - 1], nval[j - 1])});
        minwrk = max({minwrk, 3 * (mval[j - 1] + nval[j - 1]), mval[j - 1] * (mval[j - 1] + max({mval[j - 1], nval[j - 1], nrhs}) + 1) + nval[j - 1] * min(nval[j - 1], mval[j - 1])});
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
    } else if (nrhs < 0) {
        info = -6;
    } else if (lda < mmax) {
        info = -11;
    } else if (ldx < mmax) {
        info = -17;
    } else if (ldq < mmax) {
        info = -21;
    } else if (ldpt < mnmax) {
        info = -23;
    } else if (minwrk > lwork) {
        info = -27;
    }
    //
    if (info != 0) {
        Mxerbla("Rchkbd", -info);
        return;
    }
    //
    //     Initialize constants
    //
    path[0] = 'D';
    path[1] = 'B';
    path[2] = 'D';
    path[3] = '\0';
    nfail = 0;
    ntest = 0;
    unfl = Rlamch("Safe minimum");
#if defined ___MPLAPACK_DEBUG_COMPARE_WITH_QD___
    unfl = unfl * 1e+16
#endif
    ovfl = Rlamch("Overflow");
    ulp = Rlamch("Precision");
    ulpinv = one / ulp;
    log2ui = castINTEGER(log(ulpinv) / log(two));
    rtunfl = sqrt(unfl);
    rtovfl = sqrt(ovfl);
    abstol = castREAL(2) * unfl;
    //
    //     Loop over sizes, types
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        m = mval[jsize - 1];
        n = nval[jsize - 1];
        mnmin = min(m, n);
        amninv = one / max({m, n, (INTEGER)1});
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_290;
            }
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            for (j = 1; j <= 34; j = j + 1) {
                result[j - 1] = -one;
            }
            //
            uplo = ' ';
            //
            //           Compute "A"
            //
            //           Control parameters:
            //
            //           KMAGN  KMODE        KTYPE
            //       =1  O(1)   clustered 1  zero
            //       =2  large  clustered 2  identity
            //       =3  small  exponential  (none)
            //       =4         arithmetic   diagonal, (w/ eigenvalues)
            //       =5         random       symmetric, w/ eigenvalues
            //       =6                      nonsymmetric, w/ singular values
            //       =7                      random diagonal
            //       =8                      random symmetric
            //       =9                      random nonsymmetric
            //       =10                     random bidiagonal (log. distrib.)
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
            anorm = (rtovfl * ulp) * amninv;
            goto statement_70;
        //
        statement_60:
            anorm = rtunfl * castREAL(max(m, n)) * ulpinv;
            goto statement_70;
        //
        statement_70:
            //
            Rlaset("Full", lda, n, zero, zero, a, lda);
            iinfo = 0;
            cond = ulpinv;
            //
            bidiag = false;
            if (itype == 1) {
                //
                //              Zero matrix
                //
                iinfo = 0;
                //
            } else if (itype == 2) {
                //
                //              Identity
                //
                for (jcol = 1; jcol <= mnmin; jcol = jcol + 1) {
                    a[(jcol - 1) + (jcol - 1) * lda] = anorm;
                }
                //
            } else if (itype == 4) {
                //
                //              Diagonal Matrix, [Eigen]values Specified
                //
                Rlatms(mnmin, mnmin, "S", iseed, "N", work, imode, cond, anorm, 0, 0, "N", a, lda, &work[(mnmin + 1) - 1], iinfo);
                //
            } else if (itype == 5) {
                //
                //              Symmetric, eigenvalues specified
                //
                Rlatms(mnmin, mnmin, "S", iseed, "S", work, imode, cond, anorm, m, n, "N", a, lda, &work[(mnmin + 1) - 1], iinfo);
                //
            } else if (itype == 6) {
                //
                //              Nonsymmetric, singular values specified
                //
                Rlatms(m, n, "S", iseed, "N", work, imode, cond, anorm, m, n, "N", a, lda, &work[(mnmin + 1) - 1], iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random entries
                //
                Rlatmr(mnmin, mnmin, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(2 * mnmin + 1) - 1], 1, one, "N", iwork, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Symmetric, random entries
                //
                Rlatmr(mnmin, mnmin, "S", iseed, "S", work, 6, one, one, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(m + mnmin + 1) - 1], 1, one, "N", iwork, m, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              Nonsymmetric, random entries
                //
                Rlatmr(m, n, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(m + mnmin + 1) - 1], 1, one, "N", iwork, m, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 10) {
                //
                //              Bidiagonal, random entries
                //
                temp1 = -two * log(ulp);
                for (j = 1; j <= mnmin; j = j + 1) {
                    bd[j - 1] = exp(temp1 * Rlarnd(2, iseed));
                    if (j < mnmin) {
                        be[j - 1] = exp(temp1 * Rlarnd(2, iseed));
                    }
                }
                //
                iinfo = 0;
                bidiag = true;
                if (m >= n) {
                    uplo = 'U';
                } else {
                    uplo = 'L';
                }
            } else {
                iinfo = 1;
            }
            //
            if (iinfo == 0) {
                //
                //              Generate Right-Hand Side
                //
                if (bidiag) {
                    Rlatmr(mnmin, nrhs, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(2 * mnmin + 1) - 1], 1, one, "N", iwork, mnmin, nrhs, zero, one, "NO", y, ldx, iwork, iinfo);
                } else {
                    Rlatmr(m, nrhs, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(m + 1) - 1], 1, one, &work[(2 * m + 1) - 1], 1, one, "N", iwork, m, nrhs, zero, one, "NO", x, ldx, iwork, iinfo);
                }
            }
            //
            //           Error Exit
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Generator", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                return;
            }
        //
        statement_100:
            //
            //           Call Rgebrd and Rorgbr to compute B, Q, and P, do tests.
            //
            if (!bidiag) {
                //
                //              Compute transformations to reduce A to bidiagonal form:
                //              B := Q' * A * P.
                //
                Rlacpy(" ", m, n, a, lda, q, ldq);
                Rgebrd(m, n, q, ldq, bd, be, work, &work[(mnmin + 1) - 1], &work[(2 * mnmin + 1) - 1], lwork - 2 * mnmin, iinfo);
                //
                //              Check error code from Rgebrd.
                //
                if (iinfo != 0) {
                    write(nout, format_9998), "Rgebrd", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                Rlacpy(" ", m, n, q, ldq, pt, ldpt);
                if (m >= n) {
                    uplo = 'U';
                } else {
                    uplo = 'L';
                }
                //
                //              Generate Q
                //
                mq = m;
                if (nrhs <= 0) {
                    mq = mnmin;
                }
                Rorgbr("Q", m, mq, n, q, ldq, work, &work[(2 * mnmin + 1) - 1], lwork - 2 * mnmin, iinfo);
                //
                //              Check error code from Rorgbr.
                //
                if (iinfo != 0) {
                    write(nout, format_9998), "Rorgbr(Q)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                //              Generate P'
                //
                Rorgbr("P", mnmin, n, m, pt, ldpt, &work[(mnmin + 1) - 1], &work[(2 * mnmin + 1) - 1], lwork - 2 * mnmin, iinfo);
                //
                //              Check error code from Rorgbr.
                //
                if (iinfo != 0) {
                    write(nout, format_9998), "Rorgbr(P)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                //              Apply Q' to an M by NRHS matrix X:  Y := Q' * X.
                //
                Rgemm("Transpose", "No transpose", m, nrhs, m, one, q, ldq, x, ldx, zero, y, ldx);
                //
                //              Test 1:  Check the decomposition A := Q * B * PT
                //                   2:  Check the orthogonality of Q
                //                   3:  Check the orthogonality of PT
                //
                Rbdt01(m, n, 1, a, lda, q, ldq, bd, be, pt, ldpt, work, result[1 - 1]);
                Rort01("Columns", m, mq, q, ldq, work, lwork, result[2 - 1]);
                Rort01("Rows", mnmin, n, pt, ldpt, work, lwork, result[3 - 1]);
            }
            //
            //           Use Rbdsqr to form the SVD of the bidiagonal matrix B:
            //           B := U * S1 * VT, and compute Z = U' * Y.
            //
            Rcopy(mnmin, bd, 1, s1, 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, work, 1);
            }
            Rlacpy(" ", m, nrhs, y, ldx, z, ldx);
            Rlaset("Full", mnmin, mnmin, zero, one, u, ldpt);
            Rlaset("Full", mnmin, mnmin, zero, one, vt, ldpt);
            //
            Rbdsqr(&uplo, mnmin, mnmin, mnmin, nrhs, s1, work, vt, ldpt, u, ldpt, z, ldx, &work[(mnmin + 1) - 1], iinfo);
            //
            //           Check error code from Rbdsqr.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsqr(vects)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[4 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            //           Use Rbdsqr to compute only the singular values of the
            //           bidiagonal matrix B;  U, VT, and Z should not be modified.
            //
            Rcopy(mnmin, bd, 1, s2, 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, work, 1);
            }
            //
            Rbdsqr(&uplo, mnmin, 0, 0, 0, s2, work, vt, ldpt, u, ldpt, z, ldx, &work[(mnmin + 1) - 1], iinfo);
            //
            //           Check error code from Rbdsqr.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsqr(values)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[9 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            //           Test 4:  Check the decomposition B := U * S1 * VT
            //                5:  Check the computation Z := U' * Y
            //                6:  Check the orthogonality of U
            //                7:  Check the orthogonality of VT
            //
            Rbdt03(&uplo, mnmin, 1, bd, be, u, ldpt, s1, vt, ldpt, work, result[4 - 1]);
            Rbdt02(mnmin, nrhs, y, ldx, z, ldx, u, ldpt, work, result[5 - 1]);
            Rort01("Columns", mnmin, mnmin, u, ldpt, work, lwork, result[6 - 1]);
            Rort01("Rows", mnmin, mnmin, vt, ldpt, work, lwork, result[7 - 1]);
            //
            //           Test 8:  Check that the singular values are sorted in
            //                    non-increasing order and are non-negative
            //
            result[8 - 1] = zero;
            for (i = 1; i <= mnmin - 1; i = i + 1) {
                if (s1[i - 1] < s1[(i + 1) - 1]) {
                    result[8 - 1] = ulpinv;
                }
                if (s1[i - 1] < zero) {
                    result[8 - 1] = ulpinv;
                }
            }
            if (mnmin >= 1) {
                if (s1[mnmin - 1] < zero) {
                    result[8 - 1] = ulpinv;
                }
            }
            //
            //           Test 9:  Compare Rbdsqr with and without singular vectors
            //
            temp2 = zero;
            //
            for (j = 1; j <= mnmin; j = j + 1) {
                temp1 = abs(s1[j - 1] - s2[j - 1]) / max(REAL(sqrt(unfl) * max(s1[1 - 1], one)), REAL(ulp * max(abs(s1[j - 1]), abs(s2[j - 1]))));
                temp2 = max(temp1, temp2);
            }
            //
            result[9 - 1] = temp2;
            //
            //           Test 10:  Sturm sequence test of singular values
            //                     Go up by factors of two until it succeeds
            //
            temp1 = thresh * (half - ulp);
            //
            for (j = 0; j <= log2ui; j = j + 1) {
                //               CALL Rsvdch( MNMIN, BD, BE, S1, TEMP1, IINFO )
                if (iinfo == 0) {
                    goto statement_140;
                }
                temp1 = temp1 * two;
            }
        //
        statement_140:
            result[10 - 1] = temp1;
            //
            //           Use Rbdsqr to form the decomposition A := (QU) S (VT PT)
            //           from the bidiagonal form A := Q B PT.
            //
            if (!bidiag) {
                Rcopy(mnmin, bd, 1, s2, 1);
                if (mnmin > 0) {
                    Rcopy(mnmin - 1, be, 1, work, 1);
                }
                //
                Rbdsqr(&uplo, mnmin, n, m, nrhs, s2, work, pt, ldpt, q, ldq, y, ldx, &work[(mnmin + 1) - 1], iinfo);
                //
                //              Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
                //                   12:  Check the computation Z := U' * Q' * X
                //                   13:  Check the orthogonality of Q*U
                //                   14:  Check the orthogonality of VT*PT
                //
                Rbdt01(m, n, 0, a, lda, q, ldq, s2, dumma, pt, ldpt, work, result[11 - 1]);
                Rbdt02(m, nrhs, x, ldx, y, ldx, q, ldq, work, result[12 - 1]);
                Rort01("Columns", m, mq, q, ldq, work, lwork, result[13 - 1]);
                Rort01("Rows", mnmin, n, pt, ldpt, work, lwork, result[14 - 1]);
            }
            //
            //           Use Rbdsdc to form the SVD of the bidiagonal matrix B:
            //           B := U * S1 * VT
            //
            Rcopy(mnmin, bd, 1, s1, 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, work, 1);
            }
            Rlaset("Full", mnmin, mnmin, zero, one, u, ldpt);
            Rlaset("Full", mnmin, mnmin, zero, one, vt, ldpt);
            //
            Rbdsdc(&uplo, "I", mnmin, s1, work, u, ldpt, vt, ldpt, dum, idum, &work[(mnmin + 1) - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsdc.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsdc(vects)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[15 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            //           Use Rbdsdc to compute only the singular values of the
            //           bidiagonal matrix B;  U and VT should not be modified.
            //
            Rcopy(mnmin, bd, 1, s2, 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, work, 1);
            }
            //
            Rbdsdc(&uplo, "N", mnmin, s2, work, dum, 1, dum, 1, dum, idum, &work[(mnmin + 1) - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsdc.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsdc(values)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[18 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            //           Test 15:  Check the decomposition B := U * S1 * VT
            //                16:  Check the orthogonality of U
            //                17:  Check the orthogonality of VT
            //
            Rbdt03(&uplo, mnmin, 1, bd, be, u, ldpt, s1, vt, ldpt, work, result[15 - 1]);
            Rort01("Columns", mnmin, mnmin, u, ldpt, work, lwork, result[16 - 1]);
            Rort01("Rows", mnmin, mnmin, vt, ldpt, work, lwork, result[17 - 1]);
            //
            //           Test 18:  Check that the singular values are sorted in
            //                     non-increasing order and are non-negative
            //
            result[18 - 1] = zero;
            for (i = 1; i <= mnmin - 1; i = i + 1) {
                if (s1[i - 1] < s1[(i + 1) - 1]) {
                    result[18 - 1] = ulpinv;
                }
                if (s1[i - 1] < zero) {
                    result[18 - 1] = ulpinv;
                }
            }
            if (mnmin >= 1) {
                if (s1[mnmin - 1] < zero) {
                    result[18 - 1] = ulpinv;
                }
            }
            //
            //           Test 19:  Compare Rbdsqr with and without singular vectors
            //
            temp2 = zero;
            //
            for (j = 1; j <= mnmin; j = j + 1) {
                temp1 = abs(s1[j - 1] - s2[j - 1]) / max(REAL(sqrt(unfl) * max(s1[1 - 1], one)), REAL(ulp * max(abs(s1[1 - 1]), abs(s2[1 - 1]))));
                temp2 = max(temp1, temp2);
            }
            //
            result[19 - 1] = temp2;
            //
            //           Use Rbdsvdx to compute the SVD of the bidiagonal matrix B:
            //           B := U * S1 * VT
            //
            if (jtype == 10 || jtype == 16) {
                //              =================================
                //              Matrix types temporarily disabled
                //              =================================
                for (int p = 19; p <= 33; p++)
                    result[p] = zero;
                goto statement_270;
            }
            //
            iwbs = 1;
            iwbd = iwbs + mnmin;
            iwbe = iwbd + mnmin;
            iwbz = iwbe + mnmin;
            iwwork = iwbz + 2 * mnmin * (mnmin + 1);
            mnmin2 = max((INTEGER)1, mnmin * 2);
            //
            Rcopy(mnmin, bd, 1, &work[iwbd - 1], 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, &work[iwbe - 1], 1);
            }
            //
            Rbdsvdx(&uplo, "V", "A", mnmin, &work[iwbd - 1], &work[iwbe - 1], zero, zero, 0, 0, ns1, s1, &work[iwbz - 1], mnmin2, &work[iwwork - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsvdx.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsvdx(vects,A)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[20 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            j = iwbz;
            for (i = 1; i <= ns1; i = i + 1) {
                Rcopy(mnmin, &work[j - 1], 1, &u[(i - 1) * ldu], 1);
                j += mnmin;
                Rcopy(mnmin, &work[j - 1], 1, &vt[(i - 1)], ldpt);
                j += mnmin;
            }
            //
            //           Use Rbdsvdx to compute only the singular values of the
            //           bidiagonal matrix B;  U and VT should not be modified.
            //
            if (jtype == 9) {
                //              =================================
                //              Matrix types temporarily disabled
                //              =================================
                result[24 - 1] = zero;
                goto statement_270;
            }
            //
            Rcopy(mnmin, bd, 1, &work[iwbd - 1], 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, &work[iwbe - 1], 1);
            }
            //
            Rbdsvdx(&uplo, "N", "A", mnmin, &work[iwbd - 1], &work[iwbe - 1], zero, zero, 0, 0, ns2, s2, &work[iwbz - 1], mnmin2, &work[iwwork - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsvdx.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsvdx(values,A)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[24 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            //           Save S1 for tests 30-34.
            //
            Rcopy(mnmin, s1, 1, &work[iwbs - 1], 1);
            //
            //           Test 20:  Check the decomposition B := U * S1 * VT
            //                21:  Check the orthogonality of U
            //                22:  Check the orthogonality of VT
            //                23:  Check that the singular values are sorted in
            //                     non-increasing order and are non-negative
            //                24:  Compare Rbdsvdx with and without singular vectors
            //
            Rbdt03(&uplo, mnmin, 1, bd, be, u, ldpt, s1, vt, ldpt, &work[(iwbs + mnmin) - 1], result[20 - 1]);
            Rort01("Columns", mnmin, mnmin, u, ldpt, &work[(iwbs + mnmin) - 1], lwork - mnmin, result[21 - 1]);
            Rort01("Rows", mnmin, mnmin, vt, ldpt, &work[(iwbs + mnmin) - 1], lwork - mnmin, result[22 - 1]);
            //
            result[23 - 1] = zero;
            for (i = 1; i <= mnmin - 1; i = i + 1) {
                if (s1[i - 1] < s1[(i + 1) - 1]) {
                    result[23 - 1] = ulpinv;
                }
                if (s1[i - 1] < zero) {
                    result[23 - 1] = ulpinv;
                }
            }
            if (mnmin >= 1) {
                if (s1[mnmin - 1] < zero) {
                    result[23 - 1] = ulpinv;
                }
            }
            //
            temp2 = zero;
            for (j = 1; j <= mnmin; j = j + 1) {
                temp1 = abs(s1[j - 1] - s2[j - 1]) / max(REAL(sqrt(unfl) * max(s1[1 - 1], one)), REAL(ulp * max(abs(s1[1 - 1]), abs(s2[1 - 1]))));
                temp2 = max(temp1, temp2);
            }
            result[24 - 1] = temp2;
            anorm = s1[1 - 1];
            //
            //           Use Rbdsvdx with RANGE='I': choose random values for IL and
            //           IU, and ask for the IL-th through IU-th singular values
            //           and corresponding vectors.
            //
            for (i = 1; i <= 4; i = i + 1) {
                iseed2[i - 1] = iseed[i - 1];
            }
            if (mnmin <= 1) {
                il = 1;
                iu = mnmin;
            } else {
                il = 1 + castINTEGER((mnmin - 1) * Rlarnd(1, iseed2));
                iu = 1 + castINTEGER((mnmin - 1) * Rlarnd(1, iseed2));
                if (iu < il) {
                    itemp = iu;
                    iu = il;
                    il = itemp;
                }
            }
            //
            Rcopy(mnmin, bd, 1, &work[iwbd - 1], 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, &work[iwbe - 1], 1);
            }
            //
            Rbdsvdx(&uplo, "V", "I", mnmin, &work[iwbd - 1], &work[iwbe - 1], zero, zero, il, iu, ns1, s1, &work[iwbz - 1], mnmin2, &work[iwwork - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsvdx.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsvdx(vects,I)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[25 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            j = iwbz;
            for (i = 1; i <= ns1; i = i + 1) {
                Rcopy(mnmin, &work[j - 1], 1, &u[(i - 1) * ldu], 1);
                j += mnmin;
                Rcopy(mnmin, &work[j - 1], 1, &vt[(i - 1)], ldpt);
                j += mnmin;
            }
            //
            //           Use Rbdsvdx to compute only the singular values of the
            //           bidiagonal matrix B;  U and VT should not be modified.
            //
            Rcopy(mnmin, bd, 1, &work[iwbd - 1], 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, &work[iwbe - 1], 1);
            }
            //
            Rbdsvdx(&uplo, "N", "I", mnmin, &work[iwbd - 1], &work[iwbe - 1], zero, zero, il, iu, ns2, s2, &work[iwbz - 1], mnmin2, &work[iwwork - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsvdx.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsvdx(values,I)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[29 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            //           Test 25:  Check S1 - U' * B * VT'
            //                26:  Check the orthogonality of U
            //                27:  Check the orthogonality of VT
            //                28:  Check that the singular values are sorted in
            //                     non-increasing order and are non-negative
            //                29:  Compare Rbdsvdx with and without singular vectors
            //
            Rbdt04(&uplo, mnmin, bd, be, s1, ns1, u, ldpt, vt, ldpt, &work[(iwbs + mnmin) - 1], result[25 - 1]);
            Rort01("Columns", mnmin, ns1, u, ldpt, &work[(iwbs + mnmin) - 1], lwork - mnmin, result[26 - 1]);
            Rort01("Rows", ns1, mnmin, vt, ldpt, &work[(iwbs + mnmin) - 1], lwork - mnmin, result[27 - 1]);
            //
            result[28 - 1] = zero;
            for (i = 1; i <= ns1 - 1; i = i + 1) {
                if (s1[i - 1] < s1[(i + 1) - 1]) {
                    result[28 - 1] = ulpinv;
                }
                if (s1[i - 1] < zero) {
                    result[28 - 1] = ulpinv;
                }
            }
            if (ns1 >= 1) {
                if (s1[ns1 - 1] < zero) {
                    result[28 - 1] = ulpinv;
                }
            }
            //
            temp2 = zero;
            for (j = 1; j <= ns1; j = j + 1) {
                temp1 = abs(s1[j - 1] - s2[j - 1]) / max(REAL(sqrt(unfl) * max(s1[1 - 1], one)), REAL(ulp * max(abs(s1[1 - 1]), abs(s2[1 - 1]))));
                temp2 = max(temp1, temp2);
            }
            result[29 - 1] = temp2;
            //
            //           Use Rbdsvdx with RANGE='V': determine the values VL and VU
            //           of the IL-th and IU-th singular values and ask for all
            //           singular values in this range.
            //
            Rcopy(mnmin, &work[iwbs - 1], 1, s1, 1);
            //
            if (mnmin > 0) {
                if (il != 1) {
                    vu = s1[il - 1] + max({REAL(half * abs(s1[il - 1] - s1[(il - 1) - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                } else {
                    vu = s1[1 - 1] + max({REAL(half * abs(s1[mnmin - 1] - s1[1 - 1])), REAL(ulp * anorm), REAL(two * rtunfl)});
                }
                if (iu != ns1) {
                    vl = s1[iu - 1] - max({REAL(ulp * anorm), REAL(two * rtunfl), REAL(half * abs(s1[(iu + 1) - 1] - s1[iu - 1]))});
                } else {
                    vl = s1[ns1 - 1] - max({REAL(ulp * anorm), REAL(two * rtunfl), REAL(half * abs(s1[mnmin - 1] - s1[1 - 1]))});
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
            //
            Rcopy(mnmin, bd, 1, &work[iwbd - 1], 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, &work[iwbe - 1], 1);
            }
            //
            Rbdsvdx(&uplo, "V", "V", mnmin, &work[iwbd - 1], &work[iwbe - 1], vl, vu, 0, 0, ns1, s1, &work[iwbz - 1], mnmin2, &work[iwwork - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsvdx.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsvdx(vects,V)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[30 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            j = iwbz;
            for (i = 1; i <= ns1; i = i + 1) {
                Rcopy(mnmin, &work[j - 1], 1, &u[(i - 1) * ldu], 1);
                j += mnmin;
                Rcopy(mnmin, &work[j - 1], 1, &vt[(i - 1)], ldpt);
                j += mnmin;
            }
            //
            //           Use Rbdsvdx to compute only the singular values of the
            //           bidiagonal matrix B;  U and VT should not be modified.
            //
            Rcopy(mnmin, bd, 1, &work[iwbd - 1], 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, &work[iwbe - 1], 1);
            }
            //
            Rbdsvdx(&uplo, "N", "V", mnmin, &work[iwbd - 1], &work[iwbe - 1], vl, vu, 0, 0, ns2, s2, &work[iwbz - 1], mnmin2, &work[iwwork - 1], iwork, iinfo);
            //
            //           Check error code from Rbdsvdx.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Rbdsvdx(values,V)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[34 - 1] = ulpinv;
                    goto statement_270;
                }
            }
            //
            //           Test 30:  Check S1 - U' * B * VT'
            //                31:  Check the orthogonality of U
            //                32:  Check the orthogonality of VT
            //                33:  Check that the singular values are sorted in
            //                     non-increasing order and are non-negative
            //                34:  Compare Rbdsvdx with and without singular vectors
            //
            Rbdt04(&uplo, mnmin, bd, be, s1, ns1, u, ldpt, vt, ldpt, &work[(iwbs + mnmin) - 1], result[30 - 1]);
            Rort01("Columns", mnmin, ns1, u, ldpt, &work[(iwbs + mnmin) - 1], lwork - mnmin, result[31 - 1]);
            Rort01("Rows", ns1, mnmin, vt, ldpt, &work[(iwbs + mnmin) - 1], lwork - mnmin, result[32 - 1]);
            //
            result[33 - 1] = zero;
            for (i = 1; i <= ns1 - 1; i = i + 1) {
                if (s1[i - 1] < s1[(i + 1) - 1]) {
                    result[28 - 1] = ulpinv;
                }
                if (s1[i - 1] < zero) {
                    result[28 - 1] = ulpinv;
                }
            }
            if (ns1 >= 1) {
                if (s1[ns1 - 1] < zero) {
                    result[28 - 1] = ulpinv;
                }
            }
            //
            temp2 = zero;
            for (j = 1; j <= ns1; j = j + 1) {
                temp1 = abs(s1[j - 1] - s2[j - 1]) / max(REAL(sqrt(unfl) * max(s1[1 - 1], one)), REAL(ulp * max(abs(s1[1 - 1]), abs(s2[1 - 1]))));
                temp2 = max(temp1, temp2);
            }
            result[34 - 1] = temp2;
        //
        //           End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_270:
            //
            for (j = 1; j <= 34; j = j + 1) {
                if (result[j - 1] >= thresh) {
                    if (nfail == 0) {
                        Rlahd2(nout, path);
                    }
                    sprintnum_short(buf, result[j - 1]);
                    write(nout, "(' M=',i5,', N=',i5,', type ',i2,', seed=',4(i4,','),' test(',i2,"
                                "')=',a)"),
                        m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3], j, buf;
                    nfail++;
                }
            }
            if (!bidiag) {
                ntest += 34;
            } else {
                ntest += 30;
            }
        //
        statement_290:;
        }
    }
    //
    //     Summary
    //
    Alasum(path, nout, nfail, ntest, 0);
    //
    //     End of Rchkbd
    //
}

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

#if defined ___MPLAPACK_DEBUG_COMPARE_WITH_DOUBLE___
#include <lapacke.h>
#endif

void Cchkbd(INTEGER const nsizes, INTEGER *mval, INTEGER *nval, INTEGER const ntypes, bool *dotype, INTEGER const nrhs, INTEGER *iseed, REAL const thresh, COMPLEX *a, INTEGER const lda, REAL *bd, REAL *be, REAL *s1, REAL *s2, COMPLEX *x, INTEGER const ldx, COMPLEX *y, COMPLEX *z, COMPLEX *q, INTEGER const ldq, COMPLEX *pt, INTEGER const ldpt, COMPLEX *u, COMPLEX *vt, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const nout, INTEGER &info) {
    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 16;
    INTEGER ktype[16] = {1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9, 10};
    INTEGER kmagn[16] = {1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 0};
    INTEGER kmode[16] = {0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0};
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
    INTEGER jsize = 0;
    INTEGER m = 0;
    INTEGER n = 0;
    INTEGER mnmin = 0;
    REAL amninv = 0.0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ioldsd[4];
    REAL result[14];
    char uplo;
    INTEGER itype = 0;
    INTEGER imode = 0;
    REAL anorm = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    bool bidiag = false;
    INTEGER jcol = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER iwork[1];
    const REAL zero = 0.0;
    REAL temp1 = 0.0;
    INTEGER mq = 0;
    INTEGER i = 0;
    REAL temp2 = 0.0;
    const REAL half = 0.5e0;
    REAL dumma[1];
    INTEGER ldvt = ldpt;
    INTEGER ldu = ldpt;
    INTEGER ldy = ldx;
    INTEGER ldz = ldx;
    char buf[1024];
    static const char *format_9998 = "(' Cchkbd: ',a,' returned INFO=',i6,'.',/,9x,'M=',i6,', N=',i6,"
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
        Mxerbla("Cchkbd", -info);
        return;
    }
    //
    //     Initialize constants
    //
    path[0] = 'C';
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
                goto statement_170;
            }
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            for (j = 1; j <= 14; j = j + 1) {
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
            Claset("Full", lda, n, czero, czero, a, lda);
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
                Clatms(mnmin, mnmin, "S", iseed, "N", rwork, imode, cond, anorm, 0, 0, "N", a, lda, work, iinfo);
                //
            } else if (itype == 5) {
                //
                //              Symmetric, eigenvalues specified
                //
                Clatms(mnmin, mnmin, "S", iseed, "S", rwork, imode, cond, anorm, m, n, "N", a, lda, work, iinfo);
                //
            } else if (itype == 6) {
                //
                //              Nonsymmetric, singular values specified
                //
                Clatms(m, n, "S", iseed, "N", rwork, imode, cond, anorm, m, n, "N", a, lda, work, iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random entries
                //
                Clatmr(mnmin, mnmin, "S", iseed, "N", work, 6, one, cone, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(2 * mnmin + 1) - 1], 1, one, "N", iwork, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Symmetric, random entries
                //
                Clatmr(mnmin, mnmin, "S", iseed, "S", work, 6, one, cone, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(m + mnmin + 1) - 1], 1, one, "N", iwork, m, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              Nonsymmetric, random entries
                //
                Clatmr(m, n, "S", iseed, "N", work, 6, one, cone, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(m + mnmin + 1) - 1], 1, one, "N", iwork, m, n, zero, anorm, "NO", a, lda, iwork, iinfo);
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
                    Clatmr(mnmin, nrhs, "S", iseed, "N", work, 6, one, cone, "T", "N", &work[(mnmin + 1) - 1], 1, one, &work[(2 * mnmin + 1) - 1], 1, one, "N", iwork, mnmin, nrhs, zero, one, "NO", y, ldx, iwork, iinfo);
                } else {
                    Clatmr(m, nrhs, "S", iseed, "N", work, 6, one, cone, "T", "N", &work[(m + 1) - 1], 1, one, &work[(2 * m + 1) - 1], 1, one, "N", iwork, m, nrhs, zero, one, "NO", x, ldx, iwork, iinfo);
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
            //           Call Cgebrd and Cungbr to compute B, Q, and P, do tests.
            //
            if (!bidiag) {
                //
                //              Compute transformations to reduce A to bidiagonal form:
                //              B := Q' * A * P.
                //
                Clacpy(" ", m, n, a, lda, q, ldq);
                //                printf("a="); printmat(m, n, a, lda); printf("\n");
                Cgebrd(m, n, q, ldq, bd, be, work, &work[(mnmin + 1) - 1], &work[(2 * mnmin + 1) - 1], lwork - 2 * mnmin, iinfo);
                //                printf("bd="); printvec(bd, mnmin); printf("\n");
                //                printf("be="); printvec(be, mnmin - 1); printf("\n");
#ifdef ___MPLAPACK_DEBUG_COMPARE_WITH_DOUBLE___
                {
                    __complex__ double *a_d = new __complex__ double[m * n];
                    double *bd_d = new double[mnmin];
                    double *be_d = new double[mnmin];
                    __complex__ double *tauq_d = new __complex__ double[mnmin];
                    __complex__ double *taup_d = new __complex__ double[mnmin];
                    int lda_d = m;
                    for (int pp = 0; pp < m; pp++) {
                        for (int qq = 0; qq < n; qq++) {
                            __real__ a_d[pp + qq * lda_d] = cast2double(a[pp + qq * lda].real());
                            __imag__ a_d[pp + qq * lda_d] = cast2double(a[pp + qq * lda].imag());
                        }
                    }
                    LAPACKE_zgebrd(LAPACK_COL_MAJOR, (int)m, (int)n, a_d, lda_d, bd_d, be_d, tauq_d, taup_d);
                    printf("bd_d=");
                    printvec(bd_d, mnmin);
                    printf("\n");
                    printf("be_d=");
                    printvec(be_d, mnmin - 1);
                    printf("\n");
                    delete[] taup_d;
                    delete[] tauq_d;
                    delete[] bd_d;
                    delete[] be_d;
                    delete[] a_d;
                }
#endif
                //
                //              Check error code from Cgebrd.
                //
                if (iinfo != 0) {
                    write(nout, format_9998), "Cgebrd", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                Clacpy(" ", m, n, q, ldq, pt, ldpt);
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
                Cungbr("Q", m, mq, n, q, ldq, work, &work[(2 * mnmin + 1) - 1], lwork - 2 * mnmin, iinfo);
                //
                //              Check error code from Cungbr.
                //
                if (iinfo != 0) {
                    write(nout, format_9998), "Cungbr(Q)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                //              Generate P'
                //
                Cungbr("P", mnmin, n, m, pt, ldpt, &work[(mnmin + 1) - 1], &work[(2 * mnmin + 1) - 1], lwork - 2 * mnmin, iinfo);
                //
                //              Check error code from Cungbr.
                //
                if (iinfo != 0) {
                    write(nout, format_9998), "Cungbr(P)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                    info = abs(iinfo);
                    return;
                }
                //
                //              Apply Q' to an M by NRHS matrix X:  Y := Q' * X.
                //
                Cgemm("Conjugate transpose", "No transpose", m, nrhs, m, cone, q, ldq, x, ldx, czero, y, ldx);
                //
                //              Test 1:  Check the decomposition A := Q * B * PT
                //                   2:  Check the orthogonality of Q
                //                   3:  Check the orthogonality of PT
                //
                Cbdt01(m, n, 1, a, lda, q, ldq, bd, be, pt, ldpt, work, rwork, result[1 - 1]);
                Cunt01("Columns", m, mq, q, ldq, work, lwork, rwork, result[2 - 1]);
                Cunt01("Rows", mnmin, n, pt, ldpt, work, lwork, rwork, result[3 - 1]);
            }
            //
            //           Use Cbdsqr to form the SVD of the bidiagonal matrix B:
            //           B := U * S1 * VT, and compute Z = U' * Y.
            //
            Rcopy(mnmin, bd, 1, s1, 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, rwork, 1);
            }
            Clacpy(" ", m, nrhs, y, ldx, z, ldx);
            Claset("Full", mnmin, mnmin, czero, cone, u, ldpt);
            Claset("Full", mnmin, mnmin, czero, cone, vt, ldpt);
            //
            Cbdsqr(&uplo, mnmin, mnmin, mnmin, nrhs, s1, rwork, vt, ldpt, u, ldpt, z, ldx, &rwork[(mnmin + 1) - 1], iinfo);
            //
            //           Check error code from Cbdsqr.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Cbdsqr(vects)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[4 - 1] = ulpinv;
                    goto statement_150;
                }
            }
            //
            //           Use Cbdsqr to compute only the singular values of the
            //           bidiagonal matrix B;  U, VT, and Z should not be modified.
            //
            Rcopy(mnmin, bd, 1, s2, 1);
            if (mnmin > 0) {
                Rcopy(mnmin - 1, be, 1, rwork, 1);
            }
            //
            Cbdsqr(&uplo, mnmin, 0, 0, 0, s2, rwork, vt, ldpt, u, ldpt, z, ldx, &rwork[(mnmin + 1) - 1], iinfo);
            //
            //           Check error code from Cbdsqr.
            //
            if (iinfo != 0) {
                write(nout, format_9998), "Cbdsqr(values)", iinfo, m, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(iinfo);
                if (iinfo < 0) {
                    return;
                } else {
                    result[9 - 1] = ulpinv;
                    goto statement_150;
                }
            }
            //
            //           Test 4:  Check the decomposition B := U * S1 * VT
            //                5:  Check the computation Z := U' * Y
            //                6:  Check the orthogonality of U
            //                7:  Check the orthogonality of VT
            //
            Cbdt03(&uplo, mnmin, 1, bd, be, u, ldpt, s1, vt, ldpt, work, result[4 - 1]);
            Cbdt02(mnmin, nrhs, y, ldx, z, ldx, u, ldpt, work, rwork, result[5 - 1]);
            Cunt01("Columns", mnmin, mnmin, u, ldpt, work, lwork, rwork, result[6 - 1]);
            Cunt01("Rows", mnmin, mnmin, vt, ldpt, work, lwork, rwork, result[7 - 1]);
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
            //           Test 9:  Compare Cbdsqr with and without singular vectors
            //
            temp2 = zero;
            //
            for (j = 1; j <= mnmin; j = j + 1) {
                temp1 = abs(s1[j - 1] - s2[j - 1]) / max({REAL(sqrt(unfl) * max(s1[1 - 1], one)), REAL(ulp * max(abs(s1[j - 1]), abs(s2[j - 1])))});
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
                Rsvdch(mnmin, bd, be, s1, temp1, iinfo);
                if (iinfo == 0) {
                    goto statement_140;
                }
                temp1 = temp1 * two;
            }
        //
        statement_140:
            result[10 - 1] = temp1;
            //
            //           Use Cbdsqr to form the decomposition A := (QU) S (VT PT)
            //           from the bidiagonal form A := Q B PT.
            //
            if (!bidiag) {
                Rcopy(mnmin, bd, 1, s2, 1);
                if (mnmin > 0) {
                    Rcopy(mnmin - 1, be, 1, rwork, 1);
                }
                //
                Cbdsqr(&uplo, mnmin, n, m, nrhs, s2, rwork, pt, ldpt, q, ldq, y, ldx, &rwork[(mnmin + 1) - 1], iinfo);
                //
                //              Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
                //                   12:  Check the computation Z := U' * Q' * X
                //                   13:  Check the orthogonality of Q*U
                //                   14:  Check the orthogonality of VT*PT
                //
                Cbdt01(m, n, 0, a, lda, q, ldq, s2, dumma, pt, ldpt, work, rwork, result[11 - 1]);
                Cbdt02(m, nrhs, x, ldx, y, ldx, q, ldq, work, rwork, result[12 - 1]);
                Cunt01("Columns", m, mq, q, ldq, work, lwork, rwork, result[13 - 1]);
                Cunt01("Rows", mnmin, n, pt, ldpt, work, lwork, rwork, result[14 - 1]);
            }
        //
        //           End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_150:
            for (j = 1; j <= 14; j = j + 1) {
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
                ntest += 14;
            } else {
                ntest += 5;
            }
        //
        statement_170:;
        }
    }
    //
    //     Summary
    //
    Alasum(path, nout, nfail, ntest, 0);
    //
    //     End of Cchkbd
    //
}

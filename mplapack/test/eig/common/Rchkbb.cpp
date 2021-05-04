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

void Rchkbb(INTEGER const nsizes, INTEGER *mval, INTEGER *nval, INTEGER const nwdths, INTEGER *kk, INTEGER const ntypes, bool *dotype, INTEGER const nrhs, INTEGER *iseed, REAL const thresh, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *ab, INTEGER const ldab, REAL *bd, REAL *be, REAL *q, INTEGER const ldq, REAL *p, INTEGER const ldp, REAL *c, INTEGER const ldc, REAL *cc, REAL *work, INTEGER const lwork, REAL *result, INTEGER &info) {
    FEM_CMN_SVE(Rchkbb);
    common_write write(cmn);
    const INTEGER maxtyp = 15;
    if (is_called_first_time) {
        data((values, 1, 2, 5 * datum(4), 5 * datum(6), 3 * datum(9))), ktype;
        {
            data_values data;
            data.values, 2 * datum(1), 3 * datum(1), 2, 3, 3 * datum(1), 2, 3, 1;
            data.values, 2, 3;
            data, kmagn;
        }
        {
            data_values data;
            data.values, 2 * datum(0), 4, 3, 1, 4, 4, 4, 3;
            data.values, 1, 4, 4, 0, 0, 0;
            data, kmode;
        }
    }
    INTEGER ntestt = 0;
    bool badmm = false;
    bool badnn = false;
    INTEGER mmax = 0;
    INTEGER nmax = 0;
    INTEGER mnmax = 0;
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
    INTEGER m = 0;
    INTEGER n = 0;
    INTEGER mnmin = 0;
    REAL amninv = 0.0;
    INTEGER jwidth = 0;
    INTEGER k = 0;
    INTEGER kl = 0;
    INTEGER ku = 0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ntest = 0;
    arr_1d<4, int> ioldsd;
    INTEGER itype = 0;
    INTEGER imode = 0;
    REAL anorm = 0.0;
    const REAL zero = 0.0;
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    INTEGER jcol = 0;
    arr_1d<1, int> idumma;
    INTEGER i = 0;
    INTEGER jr = 0;
    static const char *format_9999 = "(' Rchkbb: ',a,' returned INFO=',i5,'.',/,9x,'M=',i5,' N=',i5,' K=',i5,"
                                     "', JTYPE=',i5,', ISEED=(',3(i5,','),i5,')')";
    //
    //  -- LAPACK test routine (input) --
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
    badmm = false;
    badnn = false;
    mmax = 1;
    nmax = 1;
    mnmax = 1;
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
    }
    //
    badnnb = false;
    kmax = 0;
    for (j = 1; j <= nwdths; j = j + 1) {
        kmax = max(kmax, kk[j - 1]);
        if (kk[j - 1] < 0) {
            badnnb = true;
        }
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
    } else if (nwdths < 0) {
        info = -4;
    } else if (badnnb) {
        info = -5;
    } else if (ntypes < 0) {
        info = -6;
    } else if (nrhs < 0) {
        info = -8;
    } else if (lda < nmax) {
        info = -13;
    } else if (ldab < 2 * kmax + 1) {
        info = -15;
    } else if (ldq < nmax) {
        info = -19;
    } else if (ldp < nmax) {
        info = -21;
    } else if (ldc < nmax) {
        info = -23;
    } else if ((max(lda, nmax) + 1) * nmax > lwork) {
        info = -26;
    }
    //
    if (info != 0) {
        Mxerbla("Rchkbb", -info);
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
    //     Loop over sizes, widths, types
    //
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        m = mval[jsize - 1];
        n = nval[jsize - 1];
        mnmin = min(m, n);
        amninv = one / (max({(INTEGER)1, m, n})).real();
        //
        for (jwidth = 1; jwidth <= nwdths; jwidth = jwidth + 1) {
            k = kk[jwidth - 1];
            if (k >= m && k >= n) {
                goto statement_150;
            }
            kl = max({(INTEGER)0, min(m - 1, k)});
            ku = max({(INTEGER)0, min(n - 1, k)});
            //
            if (nsizes != 1) {
                mtypes = min(maxtyp, ntypes);
            } else {
                mtypes = min(maxtyp + 1, ntypes);
            }
            //
            for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
                if (!dotype[jtype - 1]) {
                    goto statement_140;
                }
                nmats++;
                ntest = 0;
                //
                for (j = 1; j <= 4; j = j + 1) {
                    ioldsd[j - 1] = iseed[j - 1];
                }
                //
                //              Compute "A".
                //
                //              Control parameters:
                //
                //                  KMAGN  KMODE        KTYPE
                //              =1  O(1)   clustered 1  zero
                //              =2  large  clustered 2  identity
                //              =3  small  exponential  (none)
                //              =4         arithmetic   diagonal, (w/ singular values)
                //              =5         random log   (none)
                //              =6         random       nonhermitian, w/ singular values
                //              =7                      (none)
                //              =8                      (none)
                //              =9                      random nonhermitian
                //
                if (mtypes > maxtyp) {
                    goto statement_90;
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
                anorm = (rtovfl * ulp) * amninv;
                goto statement_70;
            //
            statement_60:
                anorm = rtunfl * max(m, n) * ulpinv;
                goto statement_70;
            //
            statement_70:
                //
                Rlaset("Full", lda, n, zero, zero, a, lda);
                Rlaset("Full", ldab, n, zero, zero, ab, ldab);
                iinfo = 0;
                cond = ulpinv;
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
                        a[(jcol - 1) + (jcol - 1) * lda] = anorm;
                    }
                    //
                } else if (itype == 4) {
                    //
                    //                 Diagonal Matrix, singular values specified
                    //
                    dlatms(m, n, "S", iseed, "N", work, imode, cond, anorm, 0, 0, "N", a, lda, &work[(m + 1) - 1], iinfo);
                    //
                } else if (itype == 6) {
                    //
                    //                 Nonhermitian, singular values specified
                    //
                    dlatms(m, n, "S", iseed, "N", work, imode, cond, anorm, kl, ku, "N", a, lda, &work[(m + 1) - 1], iinfo);
                    //
                } else if (itype == 9) {
                    //
                    //                 Nonhermitian, random entries
                    //
                    dlatmr(m, n, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, kl, ku, zero, anorm, "N", a, lda, idumma, iinfo);
                    //
                } else {
                    //
                    iinfo = 1;
                }
                //
                //              Generate Right-Hand Side
                //
                dlatmr(m, nrhs, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(m + 1) - 1], 1, one, &work[(2 * m + 1) - 1], 1, one, "N", idumma, m, nrhs, zero, one, "NO", c, ldc, idumma, iinfo);
                //
                if (iinfo != 0) {
                    write(nounit, format_9999), "Generator", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    return;
                }
            //
            statement_90:
                //
                //              Copy A to band storage.
                //
                for (j = 1; j <= n; j = j + 1) {
                    for (i = max((INTEGER)1, j - ku); i <= min(m, j + kl); i = i + 1) {
                        ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab] = a[(i - 1) + (j - 1) * lda];
                    }
                }
                //
                //              Copy C
                //
                Rlacpy("Full", m, nrhs, c, ldc, cc, ldc);
                //
                //              Call Rgbbrd to compute B, Q and P, and to update C.
                //
                Rgbbrd("B", m, n, nrhs, kl, ku, ab, ldab, bd, be, q, ldq, p, ldp, cc, ldc, work, iinfo);
                //
                if (iinfo != 0) {
                    write(nounit, format_9999), "Rgbbrd", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    if (iinfo < 0) {
                        return;
                    } else {
                        result[1 - 1] = ulpinv;
                        goto statement_120;
                    }
                }
                //
                //              Test 1:  Check the decomposition A := Q * B * P'
                //                   2:  Check the orthogonality of Q
                //                   3:  Check the orthogonality of P
                //                   4:  Check the computation of Q' * C
                //
                Rbdt01(m, n, -1, a, lda, q, ldq, bd, be, p, ldp, work, result[1 - 1]);
                Rort01("Columns", m, m, q, ldq, work, lwork, result[2 - 1]);
                Rort01("Rows", n, n, p, ldp, work, lwork, result[3 - 1]);
                Rbdt02(m, nrhs, c, ldc, cc, ldc, q, ldq, work, result[4 - 1]);
                //
                //              End of Loop -- Check for RESULT(j) > THRESH
                //
                ntest = 4;
            statement_120:
                ntestt += ntest;
                //
                //              Print out tests which fail.
                //
                for (jr = 1; jr <= ntest; jr = jr + 1) {
                    if (result[jr - 1] >= thresh) {
                        if (nerrs == 0) {
                            Rlahd2(nounit, "DBB");
                        }
                        nerrs++;
                        write(nounit, "(' M =',i4,' N=',i4,', K=',i3,', seed=',4(i4,','),' type ',i2,"
                                      "', test(',i2,')=',g10.3)"),
                            m, n, k, ioldsd, jtype, jr, result(jr);
                    }
                }
            //
            statement_140:;
            }
        statement_150:;
        }
    }
    //
    //     Summary
    //
    Rlasum("DBB", nounit, nerrs, ntestt);
    //
    //     End of Rchkbb
    //
}

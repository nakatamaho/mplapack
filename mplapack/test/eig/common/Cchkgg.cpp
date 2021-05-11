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

void Cchkgg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, bool const tstdif, REAL const thrshn, INTEGER const nounit, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *h, COMPLEX *t, COMPLEX *s1, COMPLEX *s2, COMPLEX *p1, COMPLEX *p2, COMPLEX *u, INTEGER const ldu, COMPLEX *v, COMPLEX *q, COMPLEX *z, COMPLEX *alpha1, COMPLEX *beta1, COMPLEX *alpha3, COMPLEX *beta3, COMPLEX *evectl, COMPLEX *evectr, COMPLEX *work, INTEGER const lwork, REAL *rwork, bool *llwork, REAL *result, INTEGER &info) {
    FEM_CMN_SVE(Cchkgg);
    iseed([4]);
    a([lda * star]);
    b([lda * star]);
    h([lda * star]);
    t([lda * star]);
    s1([lda * star]);
    s2([lda * star]);
    p1([lda * star]);
    p2([lda * star]);
    u([ldu * star]);
    v([ldu * star]);
    q([ldu * star]);
    z([ldu * star]);
    evectl([ldu * star]);
    evectr([ldu * star]);
    result([15]);
    common_write write(cmn);
    INTEGER *kadd(sve.kadd, [6]);
    const INTEGER maxtyp = 26;
    INTEGER *kamagn(sve.kamagn, [maxtyp]);
    INTEGER *katype(sve.katype, [maxtyp]);
    INTEGER *kazero(sve.kazero, [maxtyp]);
    INTEGER *kbmagn(sve.kbmagn, [maxtyp]);
    INTEGER *kbtype(sve.kbtype, [maxtyp]);
    INTEGER *kbzero(sve.kbzero, [maxtyp]);
    INTEGER *kclass(sve.kclass, [maxtyp]);
    INTEGER *ktrian(sve.ktrian, [maxtyp]);
    INTEGER *kz1(sve.kz1, [6]);
    INTEGER *kz2(sve.kz2, [6]);
    bool *lasign(sve.lasign, [maxtyp]);
    bool *lbsign(sve.lbsign, [maxtyp]);
    if (is_called_first_time) {
        data((values, 15 * datum(1), 10 * datum(2), 1 * datum(3))), kclass;
        {
            static const INTEGER values[] = {0, 1, 2, 1, 3, 3};
            data_of_type<int>(FEM_VALUES_AND_SIZE), kz1;
        }
        {
            static const INTEGER values[] = {0, 0, 1, 2, 1, 1};
            data_of_type<int>(FEM_VALUES_AND_SIZE), kz2;
        }
        {
            static const INTEGER values[] = {0, 0, 0, 0, 3, 2};
            data_of_type<int>(FEM_VALUES_AND_SIZE), kadd;
        }
        {
            data_values data;
            data.values, 0, 1, 0, 1, 2, 3, 4, 1;
            data.values, 4, 4, 1, 1, 4, 4, 4, 2;
            data.values, 4, 5, 8, 7, 9, 4 * datum(4), 0;
            data, katype;
        }
        {
            data_values data;
            data.values, 0, 0, 1, 1, 2, -3, 1, 4;
            data.values, 1, 1, 4, 4, 1, 1, -4, 2;
            data.values, -4, 8 * datum(8), 0;
            data, kbtype;
        }
        {
            data_values data;
            data.values, 6 * datum(1), 2, 1, 2 * datum(2), 2 * datum(1), 2 * datum(2), 3, 1;
            data.values, 3, 4 * datum(5), 4 * datum(3), 1;
            data, kazero;
        }
        {
            data_values data;
            data.values, 6 * datum(1), 1, 2, 2 * datum(1), 2 * datum(2), 2 * datum(1), 4, 1;
            data.values, 4, 4 * datum(6), 4 * datum(4), 1;
            data, kbzero;
        }
        {
            data_values data;
            data.values, 8 * datum(1), 2, 3, 2, 3, 2, 3, 7 * datum(1);
            data.values, 2, 3, 3, 2, 1;
            data, kamagn;
        }
        {
            data_values data;
            data.values, 8 * datum(1), 3, 2, 3, 2, 2, 3, 7 * datum(1);
            data.values, 3, 2, 3, 2, 1;
            data, kbmagn;
        }
        data((values, 16 * datum(0), 10 * datum(1))), ktrian;
        {
            data_values data;
            data.values, 6 * datum(false), true, false, 2 * datum(true), 2 * datum(false), 3 * datum(true), false, true;
            data.values, 3 * datum(false), 5 * datum(true), false;
            data, lasign;
        }
        {
            data_values data;
            data.values, 7 * datum(false), true, 2 * datum(false), 2 * datum(true), 2 * datum(false), true, false, true;
            data.values, 9 * datum(false);
            data, lbsign;
        }
    }
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    INTEGER lwkopt = 0;
    const REAL zero = 0.0;
    REAL safmin = 0.0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    REAL ulpinv = 0.0;
    REAL rmagn dim1(0, 3);
    INTEGER ntestt = 0;
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER n1 = 0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ntest = 0;
    INTEGER ioldsd[4];
    INTEGER iinfo = 0;
    INTEGER in = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER iadd = 0;
    INTEGER jc = 0;
    INTEGER jr = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    COMPLEX ctemp = 0.0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER i1 = 0;
    COMPLEX cdumma[4];
    REAL dumma[4];
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    static const char *format_9998 = "(' Cchkgg: ',a,' Eigenvectors from ',a,' incorrectly ','normalized.',/,"
                                     "' Bits of error=',0p,a,',',9x,'N=',i6,', JTYPE=',i6,', ISEED=(',3(i5,"
                                     "','),i5,')')";
    static const char *format_9999 = "(' Cchkgg: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
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
    info = 0;
    //
    badnn = false;
    nmax = 1;
    for (j = 1; j <= nsizes; j = j + 1) {
        nmax = max(nmax, nn[j - 1]);
        if (nn[j - 1] < 0) {
            badnn = true;
        }
    }
    //
    lwkopt = max({(INTEGER)2 * nmax * nmax, 4 * nmax, 1});
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
        info = -10;
    } else if (ldu <= 1 || ldu < nmax) {
        info = -19;
    } else if (lwkopt > lwork) {
        info = -30;
    }
    //
    if (info != 0) {
        Mxerbla("Cchkgg", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (nsizes == 0 || ntypes == 0) {
        return;
    }
    //
    safmin = Rlamch("Safe minimum");
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    safmin = safmin / ulp;
    safmax = one / safmin;
    Rlabad(safmin, safmax);
    ulpinv = one / ulp;
    //
    //     The values RMAGN(2:3) depend on N, see below.
    //
    rmagn[0 - 1] = zero;
    rmagn[1 - 1] = one;
    //
    //     Loop over sizes, types
    //
    ntestt = 0;
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        n1 = max((INTEGER)1, n);
        rmagn[2 - 1] = safmax * ulp / n1.real();
        rmagn[3 - 1] = safmin * ulpinv * n1;
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_230;
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
            for (j = 1; j <= 15; j = j + 1) {
                result[j - 1] = zero;
            }
            //
            //           Compute A and B
            //
            //           Description of control parameters:
            //
            //           KZLASS: =1 means w/o rotation, =2 means w/ rotation,
            //                   =3 means random.
            //           KATYPE: the "type" to be passed to Clatm4 for computing A.
            //           KAZERO: the pattern of zeros on the diagonal for A:
            //                   =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
            //                   =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
            //                   =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
            //                   non-zero entries.)
            //           KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
            //                   =2: large, =3: small.
            //           LASIGN: .TRUE. if the diagonal elements of A are to be
            //                   multiplied by a random magnitude 1 number.
            //           KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
            //           KTRIAN: =0: don't fill in the upper triangle, =1: do.
            //           KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            //           RMAGN:  used to implement KAMAGN and KBMAGN.
            //
            if (mtypes > maxtyp) {
                goto statement_110;
            }
            iinfo = 0;
            if (kclass[jtype - 1] < 3) {
                //
                //              Generate A (w/o rotation)
                //
                if (abs(katype[jtype - 1]) == 3) {
                    in = 2 * ((n - 1) / 2) + 1;
                    if (in != n) {
                        Claset("Full", n, n, czero, czero, a, lda);
                    }
                } else {
                    in = n;
                }
                Clatm4(katype[jtype - 1], in, kz1[kazero[jtype - 1] - 1], kz2[kazero[jtype - 1] - 1], lasign[jtype - 1], rmagn[kamagn[jtype - 1] - 1], ulp, rmagn[(ktrian[jtype - 1] * kamagn[jtype - 1]) - 1], 4, iseed, a, lda);
                iadd = kadd[kazero[jtype - 1] - 1];
                if (iadd > 0 && iadd <= n) {
                    a[(iadd - 1) + (iadd - 1) * lda] = rmagn[kamagn[jtype - 1] - 1];
                }
                //
                //              Generate B (w/o rotation)
                //
                if (abs(kbtype[jtype - 1]) == 3) {
                    in = 2 * ((n - 1) / 2) + 1;
                    if (in != n) {
                        Claset("Full", n, n, czero, czero, b, lda);
                    }
                } else {
                    in = n;
                }
                Clatm4(kbtype[jtype - 1], in, kz1[kbzero[jtype - 1] - 1], kz2[kbzero[jtype - 1] - 1], lbsign[jtype - 1], rmagn[kbmagn[jtype - 1] - 1], one, rmagn[(ktrian[jtype - 1] * kbmagn[jtype - 1]) - 1], 4, iseed, b, lda);
                iadd = kadd[kbzero[jtype - 1] - 1];
                if (iadd != 0) {
                    b[(iadd - 1) + (iadd - 1) * ldb] = rmagn[kbmagn[jtype - 1] - 1];
                }
                //
                if (kclass[jtype - 1] == 2 && n > 0) {
                    //
                    //                 Include rotations
                    //
                    //                 Generate U, V as Householder transformations times a
                    //                 diagonal matrix.  (Note that Clarfg makes U(j,j) and
                    //                 V(j,j) real.)
                    //
                    for (jc = 1; jc <= n - 1; jc = jc + 1) {
                        for (jr = jc; jr <= n; jr = jr + 1) {
                            u[(jr - 1) + (jc - 1) * ldu] = zlarnd(3, iseed);
                            v[(jr - 1) + (jc - 1) * ldv] = zlarnd(3, iseed);
                        }
                        Clarfg(n + 1 - jc, &u[(jc - 1) + (jc - 1) * ldu], &u[((jc + 1) - 1) + (jc - 1) * ldu], 1, &work[jc - 1]);
                        work[(2 * n + jc) - 1] = sign(one, &u[(jc - 1) + (jc - 1) * ldu].real());
                        u[(jc - 1) + (jc - 1) * ldu] = cone;
                        Clarfg(n + 1 - jc, &v[(jc - 1) + (jc - 1) * ldv], &v[((jc + 1) - 1) + (jc - 1) * ldv], 1, &work[(n + jc) - 1]);
                        work[(3 * n + jc) - 1] = sign(one, &v[(jc - 1) + (jc - 1) * ldv].real());
                        v[(jc - 1) + (jc - 1) * ldv] = cone;
                    }
                    ctemp = zlarnd(3, iseed);
                    u[(n - 1) + (n - 1) * ldu] = cone;
                    work[n - 1] = czero;
                    work[(3 * n) - 1] = ctemp / abs(ctemp);
                    ctemp = zlarnd(3, iseed);
                    v[(n - 1) + (n - 1) * ldv] = cone;
                    work[(2 * n) - 1] = czero;
                    work[(4 * n) - 1] = ctemp / abs(ctemp);
                    //
                    //                 Apply the diagonal matrices
                    //
                    for (jc = 1; jc <= n; jc = jc + 1) {
                        for (jr = 1; jr <= n; jr = jr + 1) {
                            a[(jr - 1) + (jc - 1) * lda] = work[(2 * n + jr) - 1] * conj(work[(3 * n + jc) - 1]) * a[(jr - 1) + (jc - 1) * lda];
                            b[(jr - 1) + (jc - 1) * ldb] = work[(2 * n + jr) - 1] * conj(work[(3 * n + jc) - 1]) * b[(jr - 1) + (jc - 1) * ldb];
                        }
                    }
                    Cunm2r("L", "N", n, n, n - 1, u, ldu, work, a, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Cunm2r("R", "C", n, n, n - 1, v, ldu, &work[(n + 1) - 1], a, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Cunm2r("L", "N", n, n, n - 1, u, ldu, work, b, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Cunm2r("R", "C", n, n, n - 1, v, ldu, &work[(n + 1) - 1], b, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                }
            } else {
                //
                //              Random matrices
                //
                for (jc = 1; jc <= n; jc = jc + 1) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        a[(jr - 1) + (jc - 1) * lda] = rmagn[kamagn[jtype - 1] - 1] * zlarnd(4, iseed);
                        b[(jr - 1) + (jc - 1) * ldb] = rmagn[kbmagn[jtype - 1] - 1] * zlarnd(4, iseed);
                    }
                }
            }
            //
            anorm = Clange("1", n, n, a, lda, rwork);
            bnorm = Clange("1", n, n, b, lda, rwork);
        //
        statement_100:
            //
            if (iinfo != 0) {
                write(nounit, format_9999), "Generator", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                return;
            }
        //
        statement_110:
            //
            //           Call Cgeqr2, Cunm2r, and Cgghrd to compute H, T, U, and V
            //
            Clacpy(" ", n, n, a, lda, h, lda);
            Clacpy(" ", n, n, b, lda, t, lda);
            ntest = 1;
            result[1 - 1] = ulpinv;
            //
            Cgeqr2(n, n, t, lda, work, &work[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cgeqr2", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Cunm2r("L", "C", n, n, n, t, lda, work, h, lda, &work[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cunm2r", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Claset("Full", n, n, czero, cone, u, ldu);
            Cunm2r("R", "N", n, n, n, t, lda, work, u, ldu, &work[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cunm2r", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Cgghrd("V", "I", n, 1, n, h, lda, t, lda, u, ldu, v, ldu, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Cgghrd", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            ntest = 4;
            //
            //           Do tests 1--4
            //
            Cget51(1, n, a, lda, h, lda, u, ldu, v, ldu, work, rwork, result[1 - 1]);
            Cget51(1, n, b, lda, t, lda, u, ldu, v, ldu, work, rwork, result[2 - 1]);
            Cget51(3, n, b, lda, t, lda, u, ldu, u, ldu, work, rwork, result[3 - 1]);
            Cget51(3, n, b, lda, t, lda, v, ldu, v, ldu, work, rwork, result[4 - 1]);
            //
            //           Call Chgeqz to compute S1, P1, S2, P2, Q, and Z, do tests.
            //
            //           Compute T1 and UZ
            //
            //           Eigenvalues only
            //
            Clacpy(" ", n, n, h, lda, s2, lda);
            Clacpy(" ", n, n, t, lda, p2, lda);
            ntest = 5;
            result[5 - 1] = ulpinv;
            //
            Chgeqz("E", "N", "N", n, 1, n, s2, lda, p2, lda, alpha3, beta3, q, ldu, z, ldu, work, lwork, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Chgeqz(E)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            //           Eigenvalues and Full Schur Form
            //
            Clacpy(" ", n, n, h, lda, s2, lda);
            Clacpy(" ", n, n, t, lda, p2, lda);
            //
            Chgeqz("S", "N", "N", n, 1, n, s2, lda, p2, lda, alpha1, beta1, q, ldu, z, ldu, work, lwork, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Chgeqz(S)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            //           Eigenvalues, Schur Form, and Schur Vectors
            //
            Clacpy(" ", n, n, h, lda, s1, lda);
            Clacpy(" ", n, n, t, lda, p1, lda);
            //
            Chgeqz("S", "I", "I", n, 1, n, s1, lda, p1, lda, alpha1, beta1, q, ldu, z, ldu, work, lwork, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Chgeqz(V)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            ntest = 8;
            //
            //           Do Tests 5--8
            //
            Cget51(1, n, h, lda, s1, lda, q, ldu, z, ldu, work, rwork, result[5 - 1]);
            Cget51(1, n, t, lda, p1, lda, q, ldu, z, ldu, work, rwork, result[6 - 1]);
            Cget51(3, n, t, lda, p1, lda, q, ldu, q, ldu, work, rwork, result[7 - 1]);
            Cget51(3, n, t, lda, p1, lda, z, ldu, z, ldu, work, rwork, result[8 - 1]);
            //
            //           Compute the Left and Right Eigenvectors of (S1,P1)
            //
            //           9: Compute the left eigenvector Matrix without
            //              back transforming:
            //
            ntest = 9;
            result[9 - 1] = ulpinv;
            //
            //           To test "SELECT" option, compute half of the eigenvectors
            //           in one call, and half in another
            //
            i1 = n / 2;
            for (j = 1; j <= i1; j = j + 1) {
                llwork[j - 1] = true;
            }
            for (j = i1 + 1; j <= n; j = j + 1) {
                llwork[j - 1] = false;
            }
            //
            Ctgevc("L", "S", llwork, n, s1, lda, p1, lda, evectl, ldu, cdumma, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctgevc(L,S1)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            i1 = in;
            for (j = 1; j <= i1; j = j + 1) {
                llwork[j - 1] = false;
            }
            for (j = i1 + 1; j <= n; j = j + 1) {
                llwork[j - 1] = true;
            }
            //
            Ctgevc("L", "S", llwork, n, s1, lda, p1, lda, evectl[((i1 + 1) - 1) * ldevectl], ldu, cdumma, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctgevc(L,S2)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Cget52(true, n, s1, lda, p1, lda, evectl, ldu, alpha1, beta1, work, rwork, dumma[1 - 1]);
            result[9 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thrshn) {
                write(nounit, format_9998), "Left", "Ctgevc(HOWMNY=S)", dumma(2), n, jtype, ioldsd;
            }
            //
            //           10: Compute the left eigenvector Matrix with
            //               back transforming:
            //
            ntest = 10;
            result[10 - 1] = ulpinv;
            Clacpy("F", n, n, q, ldu, evectl, ldu);
            Ctgevc("L", "B", llwork, n, s1, lda, p1, lda, evectl, ldu, cdumma, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctgevc(L,B)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Cget52(true, n, h, lda, t, lda, evectl, ldu, alpha1, beta1, work, rwork, dumma[1 - 1]);
            result[10 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thrshn) {
                write(nounit, format_9998), "Left", "Ctgevc(HOWMNY=B)", dumma(2), n, jtype, ioldsd;
            }
            //
            //           11: Compute the right eigenvector Matrix without
            //               back transforming:
            //
            ntest = 11;
            result[11 - 1] = ulpinv;
            //
            //           To test "SELECT" option, compute half of the eigenvectors
            //           in one call, and half in another
            //
            i1 = n / 2;
            for (j = 1; j <= i1; j = j + 1) {
                llwork[j - 1] = true;
            }
            for (j = i1 + 1; j <= n; j = j + 1) {
                llwork[j - 1] = false;
            }
            //
            Ctgevc("R", "S", llwork, n, s1, lda, p1, lda, cdumma, ldu, evectr, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctgevc(R,S1)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            i1 = in;
            for (j = 1; j <= i1; j = j + 1) {
                llwork[j - 1] = false;
            }
            for (j = i1 + 1; j <= n; j = j + 1) {
                llwork[j - 1] = true;
            }
            //
            Ctgevc("R", "S", llwork, n, s1, lda, p1, lda, cdumma, ldu, evectr[((i1 + 1) - 1) * ldevectr], ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctgevc(R,S2)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Cget52(false, n, s1, lda, p1, lda, evectr, ldu, alpha1, beta1, work, rwork, dumma[1 - 1]);
            result[11 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thresh) {
                write(nounit, format_9998), "Right", "Ctgevc(HOWMNY=S)", dumma(2), n, jtype, ioldsd;
            }
            //
            //           12: Compute the right eigenvector Matrix with
            //               back transforming:
            //
            ntest = 12;
            result[12 - 1] = ulpinv;
            Clacpy("F", n, n, z, ldu, evectr, ldu);
            Ctgevc("R", "B", llwork, n, s1, lda, p1, lda, cdumma, ldu, evectr, ldu, n, in, work, rwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Ctgevc(R,B)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Cget52(false, n, h, lda, t, lda, evectr, ldu, alpha1, beta1, work, rwork, dumma[1 - 1]);
            result[12 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thresh) {
                write(nounit, format_9998), "Right", "Ctgevc(HOWMNY=B)", dumma(2), n, jtype, ioldsd;
            }
            //
            //           Tests 13--15 are done only on request
            //
            if (tstdif) {
                //
                //              Do Tests 13--14
                //
                Cget51(2, n, s1, lda, s2, lda, q, ldu, z, ldu, work, rwork, result[13 - 1]);
                Cget51(2, n, p1, lda, p2, lda, q, ldu, z, ldu, work, rwork, result[14 - 1]);
                //
                //              Do Test 15
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max(temp1, abs(alpha1[j - 1] - alpha3[j - 1]));
                    temp2 = max(temp2, abs(beta1[j - 1] - beta3[j - 1]));
                }
                //
                temp1 = temp1 / max({safmin, ulp * max(temp1, anorm)});
                temp2 = temp2 / max({safmin, ulp * max(temp2, bnorm)});
                result[15 - 1] = max(temp1, temp2);
                ntest = 15;
            } else {
                result[13 - 1] = zero;
                result[14 - 1] = zero;
                result[15 - 1] = zero;
                ntest = 12;
            }
        //
        //           End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_210:
            //
            ntestt += ntest;
            //
            //           Print out tests which fail.
            //
            for (jr = 1; jr <= ntest; jr = jr + 1) {
                if (result[jr - 1] >= thresh) {
                    //
                    //                 If this is the first test to fail,
                    //                 prINTEGER a header to the data file.
                    //
                    if (nerrs == 0) {
                        write(nounit, "(1x,a3,' -- Complex Generalized eigenvalue problem')"), "ZGG";
                        //
                        //                    Matrix types
                        //
                        write(nounit, "(' Matrix types (see Cchkgg for details): ')");
                        write(nounit, "(' Special Matrices:',23x,'(J''=transposed Jordan block)',/,"
                                      "'   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ',"
                                      "'6=(diag(J'',I), diag(I,J''))',/,' Diagonal Matrices:  ( ',"
                                      "'D=diag(0,1,2,...) )',/,'   7=(D,I)   9=(large*D, small*I',"
                                      "')  11=(large*I, small*D)  13=(large*D, large*I)',/,"
                                      "'   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ',"
                                      "' 14=(small*D, small*I)',/,'  15=(D, reversed D)')");
                        write(nounit, "(' Matrices Rotated by Random ',a,' Matrices U, V:',/,"
                                      "'  16=Transposed Jordan Blocks             19=geometric ',"
                                      "'alpha, beta=0,1',/,'  17=arithm. alpha&beta             ',"
                                      "'      20=arithmetic alpha, beta=0,1',/,'  18=clustered ',"
                                      "'alpha, beta=0,1            21=random alpha, beta=0,1',/,"
                                      "' Large & Small Matrices:',/,'  22=(large, small)   ',"
                                      "'23=(small,large)    24=(small,small)    25=(large,large)',/,"
                                      "'  26=random O(1) matrices.')"),
                            "Unitary";
                        //
                        //                    Tests performed
                        //
                        {
                            write_loop wloop(cmn, nounit,
                                             "(/,' Tests performed:   (H is Hessenberg, S is Schur, B, ',"
                                             "'T, P are triangular,',/,20x,'U, V, Q, and Z are ',a,"
                                             "', l and r are the',/,20x,"
                                             "'appropriate left and right eigenvectors, resp., a is',/,20x,"
                                             "'alpha, b is beta, and ',a,' means ',a,'.)',/,"
                                             "' 1 = | A - U H V',a,"
                                             "' | / ( |A| n ulp )      2 = | B - U T V',a,"
                                             "' | / ( |B| n ulp )',/,' 3 = | I - UU',a,"
                                             "' | / ( n ulp )             4 = | I - VV',a,' | / ( n ulp )',"
                                             "/,' 5 = | H - Q S Z',a,' | / ( |H| n ulp )',6x,"
                                             "'6 = | T - Q P Z',a,' | / ( |T| n ulp )',/,' 7 = | I - QQ',a,"
                                             "' | / ( n ulp )             8 = | I - ZZ',a,' | / ( n ulp )',"
                                             "/,' 9 = max | ( b S - a P )',a,"
                                             "' l | / const.  10 = max | ( b H - a T )',a,' l | / const.',"
                                             "/,' 11= max | ( b S - a P ) r | / const.   12 = max | ( b H',"
                                             "' - a T ) r | / const.',/,1x)");
                            wloop, "unitary", "*", "conjugate transpose";
                            for (j = 1; j <= 10; j = j + 1) {
                                wloop, "*";
                            }
                        }
                        //
                    }
                    nerrs++;
                    if (result[jr - 1] < 10000.0) {
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i2,' is',0p,a)"),
                            n, jtype, ioldsd, jr, result(jr);
                    } else {
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i2,' is',1p,a)"),
                            n, jtype, ioldsd, jr, result(jr);
                    }
                }
            }
        //
        statement_230:;
        }
    }
    //
    //     Summary
    //
    Rlasum("ZGG", nounit, nerrs, ntestt);
    //
    //     End of Cchkgg
    //
}

/*
 * Copyright (c) 2021-2022
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

void Cdrgev3(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *s, COMPLEX *t, COMPLEX *q, INTEGER const ldq, COMPLEX *z, COMPLEX *qe, INTEGER const ldqe, COMPLEX *alpha, COMPLEX *beta, COMPLEX *alpha1, COMPLEX *beta1, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result, INTEGER &info) {
    INTEGER ldb = lda;
    INTEGER lds = lda;
    INTEGER ldt = lda;
    INTEGER ldz = ldq;
    char buf[1024];
    common cmn;
    common_write write(cmn);
    const INTEGER maxtyp = 26;
    INTEGER kclass[26] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3};
    INTEGER kbmagn[26] = {1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1};
    INTEGER ktrian[26] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    bool lasign[26] = {false, false, false, false, false, false, true, false, true, true, false, false, true, true, true, false, true, false, false, false, true, true, true, true, true, false};
    bool lbsign[26] = {false, false, false, false, false, false, false, true, false, false, true, true, false, false, true, false, true, false, false, false, false, false, false, false, false, false};
    INTEGER kz1[6] = {0, 1, 2, 1, 3, 3};
    INTEGER kz2[6] = {0, 0, 1, 2, 1, 1};
    INTEGER kadd[6] = {0, 0, 0, 0, 3, 2};
    INTEGER katype[26] = {0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0};
    INTEGER kbtype[26] = {0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0};
    INTEGER kazero[26] = {1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1};
    INTEGER kbzero[26] = {1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1};
    INTEGER kamagn[26] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1};
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    INTEGER minwrk = 0;
    INTEGER nb = 0;
    INTEGER maxwrk = 0;
    REAL ulp = 0.0;
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    REAL ulpinv = 0.0;
    REAL rmagn[4];
    INTEGER ntestt = 0;
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER n1 = 0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ioldsd[4];
    INTEGER ierr = 0;
    INTEGER in = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER iadd = 0;
    INTEGER jc = 0;
    INTEGER jr = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    COMPLEX ctemp = 0.0;
    INTEGER i = 0;
    static const char *format_9998 = "(' Cdrgev3: ',a,' Eigenvectors from ',a,' incorrectly normalized.',/,"
                                     "' Bits of error=',0p,a,',',3x,'N=',i4,', JTYPE=',i3,', ISEED=(',3(i4,"
                                     "','),i5,')')";
    static const char *format_9999 = "(' Cdrgev3: ',a,' returned INFO=',i6,'.',/,3x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',3(i5,','),i5,')')";
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
    } else if (ldq <= 1 || ldq < nmax) {
        info = -14;
    } else if (ldqe <= 1 || ldqe < nmax) {
        info = -17;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that point in the code,
    //       as well as the preferred amount for good performance.
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv.
    //
    minwrk = 1;
    if (info == 0 && lwork >= 1) {
        minwrk = nmax * (nmax + 1);
        nb = max({(INTEGER)1, iMlaenv(1, "Cgeqrf", " ", nmax, nmax, -1, -1), iMlaenv(1, "Cunmqr", "LC", nmax, nmax, nmax, -1), iMlaenv(1, "Cungqr", " ", nmax, nmax, nmax, -1)});
        maxwrk = max({(INTEGER)2 * nmax, nmax * (nb + 1), nmax * (nmax + 1)});
        work[1 - 1] = maxwrk;
    }
    //
    if (lwork < minwrk) {
        info = -23;
    }
    //
    if (info != 0) {
        Mxerbla("Cdrgev3", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (nsizes == 0 || ntypes == 0) {
        return;
    }
    //
    ulp = Rlamch("Precision");
    safmin = Rlamch("Safe minimum");
    safmin = safmin / ulp;
    safmax = one / safmin;
    ulpinv = one / ulp;
    //
    //     The values RMAGN(2:3) depend on N, see below.
    //
    rmagn[0] = zero;
    rmagn[1] = one;
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
        rmagn[2] = safmax * ulp / castREAL(n1);
        rmagn[3] = safmin * ulpinv * castREAL(n1);
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_210;
            }
            nmats++;
            //
            //           Save ISEED in case of an error.
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            //           Generate test matrices A and B
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
            //           RMAGN: used to implement KAMAGN and KBMAGN.
            //
            if (mtypes > maxtyp) {
                goto statement_100;
            }
            ierr = 0;
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
                Clatm4(katype[jtype - 1], in, kz1[kazero[jtype - 1] - 1], kz2[kazero[jtype - 1] - 1], lasign[jtype - 1], rmagn[kamagn[jtype - 1]], ulp, rmagn[ktrian[jtype - 1] * kamagn[jtype - 1]], 2, iseed, a, lda);
                iadd = kadd[kazero[jtype - 1] - 1];
                if (iadd > 0 && iadd <= n) {
                    a[(iadd - 1) + (iadd - 1 * lda)] = rmagn[kamagn[jtype - 1]];
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
                Clatm4(kbtype[jtype - 1], in, kz1[kbzero[jtype - 1] - 1], kz2[kbzero[jtype - 1] - 1], lbsign[jtype - 1], rmagn[kbmagn[jtype - 1]], one, rmagn[ktrian[jtype - 1] * kbmagn[jtype - 1]], 2, iseed, b, lda);
                iadd = kadd[kbzero[jtype - 1] - 1];
                if (iadd != 0 && iadd <= n) {
                    b[(iadd - 1) + (iadd - 1) * ldb] = rmagn[kbmagn[jtype - 1]];
                }
                //
                if (kclass[jtype - 1] == 2 && n > 0) {
                    //
                    //                 Include rotations
                    //
                    //                 Generate Q, Z as Householder transformations times
                    //                 a diagonal matrix.
                    //
                    for (jc = 1; jc <= n - 1; jc = jc + 1) {
                        for (jr = jc; jr <= n; jr = jr + 1) {
                            q[(jr - 1) + (jc - 1) * ldq] = Clarnd(3, iseed);
                            z[(jr - 1) + (jc - 1) * ldz] = Clarnd(3, iseed);
                        }
                        Clarfg(n + 1 - jc, q[(jc - 1) + (jc - 1) * ldq], &q[((jc + 1) - 1) + (jc - 1) * ldq], 1, work[jc - 1]);
                        work[(2 * n + jc) - 1] = sign(one, q[(jc - 1) + (jc - 1) * ldq].real());
                        q[(jc - 1) + (jc - 1) * ldq] = cone;
                        Clarfg(n + 1 - jc, z[(jc - 1) + (jc - 1) * ldz], &z[((jc + 1) - 1) + (jc - 1) * ldz], 1, work[(n + jc) - 1]);
                        work[(3 * n + jc) - 1] = sign(one, z[(jc - 1) + (jc - 1) * ldz].real());
                        z[(jc - 1) + (jc - 1) * ldz] = cone;
                    }
                    ctemp = Clarnd(3, iseed);
                    q[(n - 1) + (n - 1) * ldq] = cone;
                    work[n - 1] = czero;
                    work[(3 * n) - 1] = ctemp / abs(ctemp);
                    ctemp = Clarnd(3, iseed);
                    z[(n - 1) + (n - 1) * ldz] = cone;
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
                    Cunm2r("L", "N", n, n, n - 1, q, ldq, work, a, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                    Cunm2r("R", "C", n, n, n - 1, z, ldq, &work[(n + 1) - 1], a, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                    Cunm2r("L", "N", n, n, n - 1, q, ldq, work, b, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                    Cunm2r("R", "C", n, n, n - 1, z, ldq, &work[(n + 1) - 1], b, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                }
            } else {
                //
                //              Random matrices
                //
                for (jc = 1; jc <= n; jc = jc + 1) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        a[(jr - 1) + (jc - 1) * lda] = rmagn[kamagn[jtype - 1]] * Clarnd(4, iseed);
                        b[(jr - 1) + (jc - 1) * ldb] = rmagn[kbmagn[jtype - 1]] * Clarnd(4, iseed);
                    }
                }
            }
        //
        statement_90:
            //
            if (ierr != 0) {
                write(nounit, format_9999), "Generator", ierr, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(ierr);
                return;
            }
        //
        statement_100:
            //
            for (i = 1; i <= 7; i = i + 1) {
                result[i - 1] = -one;
            }
            //
            //           Call Cggev3 to compute eigenvalues and eigenvectors.
            //
            Clacpy(" ", n, n, a, lda, s, lda);
            Clacpy(" ", n, n, b, lda, t, lda);
            Cggev3("V", "V", n, s, lda, t, lda, alpha, beta, q, ldq, z, ldq, work, lwork, rwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Cggev31", ierr, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(ierr);
                goto statement_190;
            }
            //
            //           Do the tests (1) and (2)
            //
            Cget52(true, n, a, lda, b, lda, q, ldq, alpha, beta, work, rwork, &result[1 - 1]);
            if (result[2 - 1] > thresh) {
                sprintnum_short(buf, result[2 - 1]);
                write(nounit, format_9998), "Left", "Cggev31", buf, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
            }
            //
            //           Do the tests (3) and (4)
            //
            Cget52(false, n, a, lda, b, lda, z, ldq, alpha, beta, work, rwork, &result[3 - 1]);
            if (result[4 - 1] > thresh) {
                sprintnum_short(buf, result[4 - 1]);
                write(nounit, format_9998), "Right", "Cggev31", buf, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
            }
            //
            //           Do test (5)
            //
            Clacpy(" ", n, n, a, lda, s, lda);
            Clacpy(" ", n, n, b, lda, t, lda);
            Cggev3("N", "N", n, s, lda, t, lda, alpha1, beta1, q, ldq, z, ldq, work, lwork, rwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Cggev32", ierr, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(ierr);
                goto statement_190;
            }
            //
            for (j = 1; j <= n; j = j + 1) {
                if (alpha[j - 1] != alpha1[j - 1] || beta[j - 1] != beta1[j - 1]) {
                    result[5 - 1] = ulpinv;
                }
            }
            //
            //           Do test (6): Compute eigenvalues and left eigenvectors,
            //           and test them
            //
            Clacpy(" ", n, n, a, lda, s, lda);
            Clacpy(" ", n, n, b, lda, t, lda);
            Cggev3("V", "N", n, s, lda, t, lda, alpha1, beta1, qe, ldqe, z, ldq, work, lwork, rwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Cggev33", ierr, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(ierr);
                goto statement_190;
            }
            //
            for (j = 1; j <= n; j = j + 1) {
                if (alpha[j - 1] != alpha1[j - 1] || beta[j - 1] != beta1[j - 1]) {
                    result[6 - 1] = ulpinv;
                }
            }
            //
            for (j = 1; j <= n; j = j + 1) {
                for (jc = 1; jc <= n; jc = jc + 1) {
                    if (q[(j - 1) + (jc - 1) * ldq] != qe[(j - 1) + (jc - 1) * ldqe]) {
                        result[6 - 1] = ulpinv;
                    }
                }
            }
            //
            //           Do test (7): Compute eigenvalues and right eigenvectors,
            //           and test them
            //
            Clacpy(" ", n, n, a, lda, s, lda);
            Clacpy(" ", n, n, b, lda, t, lda);
            Cggev3("N", "V", n, s, lda, t, lda, alpha1, beta1, q, ldq, qe, ldqe, work, lwork, rwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Cggev34", ierr, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
                info = abs(ierr);
                goto statement_190;
            }
            //
            for (j = 1; j <= n; j = j + 1) {
                if (alpha[j - 1] != alpha1[j - 1] || beta[j - 1] != beta1[j - 1]) {
                    result[7 - 1] = ulpinv;
                }
            }
            //
            for (j = 1; j <= n; j = j + 1) {
                for (jc = 1; jc <= n; jc = jc + 1) {
                    if (z[(j - 1) + (jc - 1) * ldz] != qe[(j - 1) + (jc - 1) * ldqe]) {
                        result[7 - 1] = ulpinv;
                    }
                }
            }
        //
        //           End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_190:
            //
            ntestt += 7;
            //
            //           Print out tests which fail.
            //
            for (jr = 1; jr <= 7; jr = jr + 1) {
                if (result[jr - 1] >= thresh) {
                    //
                    //                 If this is the first test to fail,
                    //                 print a header to the data file.
                    //
                    if (nerrs == 0) {
                        write(nounit, "(/,1x,a3,' -- Complex Generalized eigenvalue problem ',"
                                      "'driver')"),
                            "ZGV";
                        //
                        //                    Matrix types
                        //
                        write(nounit, "(' Matrix types (see Cdrgev3 for details): ')");
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
                            "Orthogonal";
                        //
                        //                    Tests performed
                        //
                        write(nounit, "(/,' Tests performed:    ',/,"
                                      "' 1 = max | ( b A - a B )''*l | / const.,',/,"
                                      "' 2 = | |VR(i)| - 1 | / ulp,',/,"
                                      "' 3 = max | ( b A - a B )*r | / const.',/,"
                                      "' 4 = | |VL(i)| - 1 | / ulp,',/,"
                                      "' 5 = 0 if W same no matter if r or l computed,',/,"
                                      "' 6 = 0 if l same no matter if l computed,',/,"
                                      "' 7 = 0 if r same no matter if r computed,',/,1x)");
                        //
                    }
                    nerrs++;
                    if (result[jr - 1] < 10000.0) {
                        sprintnum_short(buf, result[jr - 1]);
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i2,' is',0p,a)"),
                            n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3], jr, buf;
                    } else {
                        sprintnum_short(buf, result[jr - 1]);
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i2,' is',1p,a)"),
                            n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3], jr, buf;
                    }
                }
            }
        //
        statement_210:;
        }
    }
    //
    //     Summary
    //
    Alasvm("ZGV3", nounit, nerrs, ntestt, 0);
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Cdrgev3
    //
}

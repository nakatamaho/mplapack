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

void Rdrges(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *b, REAL *s, REAL *t, REAL *q, INTEGER const ldq, REAL *z, REAL *alphar, REAL *alphai, REAL *beta, REAL *work, INTEGER const lwork, REAL *result, bool *bwork, INTEGER &info) {
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
    INTEGER iasign[26] = {0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0};
    INTEGER ibsign[26] = {0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
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
    REAL safmin = 0.0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    REAL ulpinv = 0.0;
    REAL rmagn[3];
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
    INTEGER iadd = 0;
    INTEGER jc = 0;
    INTEGER jr = 0;
    INTEGER i = 0;
    INTEGER isort = 0;
    char sort;
    INTEGER rsub = 0;
    INTEGER sdim = 0;
    REAL temp1 = 0.0;
    bool ilabad = false;
    REAL temp2 = 0.0;
    INTEGER i1 = 0;
    INTEGER ierr = 0;
    INTEGER knteig = 0;
    static const char *format_9999 = "(' Rdrges: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',4(i4,','),i5,')')";
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
        minwrk = max((INTEGER)10 * (nmax + 1), 3 * nmax * nmax);
        nb = max({(INTEGER)1, iMlaenv(1, "Rgeqrf", " ", nmax, nmax, -1, -1), iMlaenv(1, "Rormqr", "LT", nmax, nmax, nmax, -1), iMlaenv(1, "Rorgqr", " ", nmax, nmax, nmax, -1)});
        maxwrk = max({(INTEGER)10 * (nmax + 1), 2 * nmax + nmax * nb, 3 * nmax * nmax});
        work[1 - 1] = maxwrk;
    }
    //
    if (lwork < minwrk) {
        info = -20;
    }
    //
    if (info != 0) {
        Mxerbla("Rdrges", -info);
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
    //     Loop over matrix sizes
    //
    ntestt = 0;
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        n1 = max((INTEGER)1, n);
        rmagn[2 - 1] = safmax * ulp / castREAL(n1);
        rmagn[3 - 1] = safmin * ulpinv * castREAL(n1);
        //
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        //        Loop over matrix types
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_180;
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
            for (j = 1; j <= 13; j = j + 1) {
                result[j - 1] = zero;
            }
            //
            //           Generate test matrices A and B
            //
            //           Description of control parameters:
            //
            //           KZLASS: =1 means w/o rotation, =2 means w/ rotation,
            //                   =3 means random.
            //           KATYPE: the "type" to be passed to Rlatm4 for computing A.
            //           KAZERO: the pattern of zeros on the diagonal for A:
            //                   =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
            //                   =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
            //                   =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
            //                   non-zero entries.)
            //           KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
            //                   =2: large, =3: small.
            //           IASIGN: 1 if the diagonal elements of A are to be
            //                   multiplied by a random magnitude 1 number, =2 if
            //                   randomly chosen diagonal blocks are to be rotated
            //                   to form 2x2 blocks.
            //           KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
            //           KTRIAN: =0: don't fill in the upper triangle, =1: do.
            //           KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            //           RMAGN: used to implement KAMAGN and KBMAGN.
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
                        Rlaset("Full", n, n, zero, zero, a, lda);
                    }
                } else {
                    in = n;
                }
                Rlatm4(katype[jtype - 1], in, kz1[kazero[jtype - 1] - 1], kz2[kazero[jtype - 1] - 1], iasign[jtype - 1], rmagn[kamagn[jtype - 1] - 1], ulp, rmagn[(ktrian[jtype - 1] * kamagn[jtype - 1]) - 1], 2, iseed, a, lda);
                iadd = kadd[kazero[jtype - 1] - 1];
                if (iadd > 0 && iadd <= n) {
                    a[(iadd - 1) + (iadd - 1) * lda] = one;
                }
                //
                //              Generate B (w/o rotation)
                //
                if (abs(kbtype[jtype - 1]) == 3) {
                    in = 2 * ((n - 1) / 2) + 1;
                    if (in != n) {
                        Rlaset("Full", n, n, zero, zero, b, lda);
                    }
                } else {
                    in = n;
                }
                Rlatm4(kbtype[jtype - 1], in, kz1[kbzero[jtype - 1] - 1], kz2[kbzero[jtype - 1] - 1], ibsign[jtype - 1], rmagn[kbmagn[jtype - 1] - 1], one, rmagn[(ktrian[jtype - 1] * kbmagn[jtype - 1]) - 1], 2, iseed, b, lda);
                iadd = kadd[kbzero[jtype - 1] - 1];
                if (iadd != 0 && iadd <= n) {
                    b[(iadd - 1) + (iadd - 1) * ldb] = one;
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
                            q[(jr - 1) + (jc - 1) * ldq] = Rlarnd(3, iseed);
                            z[(jr - 1) + (jc - 1) * ldz] = Rlarnd(3, iseed);
                        }
                        Rlarfg(n + 1 - jc, q[(jc - 1) + (jc - 1) * ldq], &q[((jc + 1) - 1) + (jc - 1) * ldq], 1, work[jc - 1]);
                        work[(2 * n + jc) - 1] = sign(one, q[(jc - 1) + (jc - 1) * ldq]);
                        q[(jc - 1) + (jc - 1) * ldq] = one;
                        Rlarfg(n + 1 - jc, z[(jc - 1) + (jc - 1) * ldz], &z[((jc + 1) - 1) + (jc - 1) * ldz], 1, work[(n + jc) - 1]);
                        work[(3 * n + jc) - 1] = sign(one, z[(jc - 1) + (jc - 1) * ldz]);
                        z[(jc - 1) + (jc - 1) * ldz] = one;
                    }
                    q[(n - 1) + (n - 1) * ldq] = one;
                    work[n - 1] = zero;
                    work[(3 * n) - 1] = sign(one, Rlarnd(2, iseed));
                    z[(n - 1) + (n - 1) * ldz] = one;
                    work[(2 * n) - 1] = zero;
                    work[(4 * n) - 1] = sign(one, Rlarnd(2, iseed));
                    //
                    //                 Apply the diagonal matrices
                    //
                    for (jc = 1; jc <= n; jc = jc + 1) {
                        for (jr = 1; jr <= n; jr = jr + 1) {
                            a[(jr - 1) + (jc - 1) * lda] = work[(2 * n + jr) - 1] * work[(3 * n + jc) - 1] * a[(jr - 1) + (jc - 1) * lda];
                            b[(jr - 1) + (jc - 1) * ldb] = work[(2 * n + jr) - 1] * work[(3 * n + jc) - 1] * b[(jr - 1) + (jc - 1) * ldb];
                        }
                    }
                    Rorm2r("L", "N", n, n, n - 1, q, ldq, work, a, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Rorm2r("R", "T", n, n, n - 1, z, ldq, &work[(n + 1) - 1], a, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Rorm2r("L", "N", n, n, n - 1, q, ldq, work, b, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Rorm2r("R", "T", n, n, n - 1, z, ldq, &work[(n + 1) - 1], b, lda, &work[(2 * n + 1) - 1], iinfo);
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
                        a[(jr - 1) + (jc - 1) * lda] = rmagn[kamagn[jtype - 1] - 1] * Rlarnd(2, iseed);
                        b[(jr - 1) + (jc - 1) * ldb] = rmagn[kbmagn[jtype - 1] - 1] * Rlarnd(2, iseed);
                    }
                }
            }
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
            for (i = 1; i <= 13; i = i + 1) {
                result[i - 1] = -one;
            }
            //
            //           Test with and without sorting of eigenvalues
            //
            for (isort = 0; isort <= 1; isort = isort + 1) {
                if (isort == 0) {
                    sort = 'N';
                    rsub = 0;
                } else {
                    sort = 'S';
                    rsub = 5;
                }
                //
                //              Call Rgges to compute H, T, Q, Z, alpha, and beta.
                //
                Rlacpy("Full", n, n, a, lda, s, lda);
                Rlacpy("Full", n, n, b, lda, t, lda);
                ntest = 1 + rsub + isort;
                result[(1 + rsub + isort) - 1] = ulpinv;
                Rgges("V", "V", &sort, Rlctes, n, s, lda, t, lda, sdim, alphar, alphai, beta, q, ldq, z, ldq, work, lwork, bwork, iinfo);
                if (iinfo != 0 && iinfo != n + 2) {
                    result[(1 + rsub + isort) - 1] = ulpinv;
                    write(nounit, format_9999), "Rgges", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    goto statement_160;
                }
                //
                ntest = 4 + rsub;
                //
                //              Do tests 1--4 (or tests 7--9 when reordering )
                //
                if (isort == 0) {
                    Rget51(1, n, a, lda, s, lda, q, ldq, z, ldq, work, result[1 - 1]);
                    Rget51(1, n, b, lda, t, lda, q, ldq, z, ldq, work, result[2 - 1]);
                } else {
                    Rget54(n, a, lda, b, lda, s, lda, t, lda, q, ldq, z, ldq, work, result[7 - 1]);
                }
                Rget51(3, n, a, lda, t, lda, q, ldq, q, ldq, work, result[(3 + rsub) - 1]);
                Rget51(3, n, b, lda, t, lda, z, ldq, z, ldq, work, result[(4 + rsub) - 1]);
                //
                //              Do test 5 and 6 (or Tests 10 and 11 when reordering):
                //              check Schur form of A and compare eigenvalues with
                //              diagonals.
                //
                ntest = 6 + rsub;
                temp1 = zero;
                //
                for (j = 1; j <= n; j = j + 1) {
                    ilabad = false;
                    if (alphai[j - 1] == zero) {
                        temp2 = (abs(alphar[j - 1] - s[(j - 1) + (j - 1) * lds]) / max({safmin, abs(alphar[j - 1]), abs(s[(j - 1) + (j - 1) * lds])}) + abs(beta[j - 1] - t[(j - 1) + (j - 1) * ldt]) / max({safmin, abs(beta[j - 1]), abs(t[(j - 1) + (j - 1) * ldt])})) / ulp;
                        //
                        if (j < n) {
                            if (s[((j + 1) - 1) + (j - 1) * lds] != zero) {
                                ilabad = true;
                                result[(5 + rsub) - 1] = ulpinv;
                            }
                        }
                        if (j > 1) {
                            if (s[(j - 1) + ((j - 1) - 1) * lds] != zero) {
                                ilabad = true;
                                result[(5 + rsub) - 1] = ulpinv;
                            }
                        }
                        //
                    } else {
                        if (alphai[j - 1] > zero) {
                            i1 = j;
                        } else {
                            i1 = j - 1;
                        }
                        if (i1 <= 0 || i1 >= n) {
                            ilabad = true;
                        } else if (i1 < n - 1) {
                            if (s[((i1 + 2) - 1) + ((i1 + 1) - 1) * lds] != zero) {
                                ilabad = true;
                                result[(5 + rsub) - 1] = ulpinv;
                            }
                        } else if (i1 > 1) {
                            if (s[(i1 - 1) + ((i1 - 1) - 1) * lds] != zero) {
                                ilabad = true;
                                result[(5 + rsub) - 1] = ulpinv;
                            }
                        }
                        if (!ilabad) {
                            Rget53(&s[(i1 - 1) + (i1 - 1) * lds], lda, &t[(i1 - 1) + (i1 - 1) * ldt], lda, beta[j - 1], alphar[j - 1], alphai[j - 1], temp2, ierr);
                            if (ierr >= 3) {
                                write(nounit, "(' Rdrges: Rget53 returned INFO=',i1,' for eigenvalue ',i6,"
                                              "'.',/,9x,'N=',i6,', JTYPE=',i6,', ISEED=(',4(i4,','),i5,"
                                              "')')"),
                                    ierr, j, n, jtype, ioldsd;
                                info = abs(ierr);
                            }
                        } else {
                            temp2 = ulpinv;
                        }
                        //
                    }
                    temp1 = max(temp1, temp2);
                    if (ilabad) {
                        write(nounit, "(' Rdrges: S not in Schur form at eigenvalue ',i6,'.',/,9x,"
                                      "'N=',i6,', JTYPE=',i6,', ISEED=(',3(i5,','),i5,')')"),
                            j, n, jtype, ioldsd;
                    }
                }
                result[(6 + rsub) - 1] = temp1;
                //
                if (isort >= 1) {
                    //
                    //                 Do test 12
                    //
                    ntest = 12;
                    result[12 - 1] = zero;
                    knteig = 0;
                    for (i = 1; i <= n; i = i + 1) {
                        if (Rlctes(alphar[i - 1], beta[i - 1]) || Rlctes(alphar[i - 1],  beta[i - 1])) {
                            knteig++;
                        }
                        if (i < n) {
			  if ((Rlctes(alphar[(i + 1) - 1], beta[(i + 1) - 1]) || Rlctes(alphar[(i + 1) - 1], beta[(i + 1) - 1])) && (!(Rlctes(alphar[i - 1], beta[i - 1]) || Rlctes(alphar[i - 1], beta[i - 1]))) && iinfo != n + 2) {
                                result[12 - 1] = ulpinv;
                            }
                        }
                    }
                    if (sdim != knteig) {
                        result[12 - 1] = ulpinv;
                    }
                }
                //
            }
        //
        //           End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_160:
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
                        write(nounit, "(/,1x,a3,' -- Real Generalized Schur form driver')"), "DGS";
                        //
                        //                    Matrix types
                        //
                        write(nounit, "(' Matrix types (see Rdrges for details): ')");
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
                        {
                            write_loop wloop(cmn, nounit,
                                             "(/,' Tests performed:  (S is Schur, T is triangular, ',"
                                             "'Q and Z are ',a,',',/,19x,"
                                             "'l and r are the appropriate left and right',/,19x,"
                                             "'eigenvectors, resp., a is alpha, b is beta, and',/,19x,a,"
                                             "' means ',a,'.)',/,' Without ordering: ',/,"
                                             "'  1 = | A - Q S Z',a,"
                                             "' | / ( |A| n ulp )      2 = | B - Q T Z',a,"
                                             "' | / ( |B| n ulp )',/,'  3 = | I - QQ',a,"
                                             "' | / ( n ulp )             4 = | I - ZZ',a,' | / ( n ulp )',"
                                             "/,'  5 = A is in Schur form S',/,"
                                             "'  6 = difference between (alpha,beta)',"
                                             "' and diagonals of (S,T)',/,' With ordering: ',/,"
                                             "'  7 = | (A,B) - Q (S,T) Z',a,' | / ( |(A,B)| n ulp )  ',/,"
                                             "'  8 = | I - QQ',a,' | / ( n ulp )            9 = | I - ZZ',"
                                             "a,' | / ( n ulp )',/,' 10 = A is in Schur form S',/,"
                                             "' 11 = difference between (alpha,beta) and diagonals',"
                                             "' of (S,T)',/,' 12 = SDIM is the correct number of ',"
                                             "'selected eigenvalues',/)");
                            wloop, "orthogonal", "'", "transpose";
                            for (j = 1; j <= 8; j = j + 1) {
                                wloop, "'";
                            }
                        }
                        //
                    }
                    nerrs++;
                    if (result[jr - 1] < 10000.0) {
		        sprintnum_short(buf, result[jr - 1]);
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i2,' is',0p,a)"),
			  n, jtype, ioldsd, jr, buf;
                    } else {
		        sprintnum_short(buf, result[jr - 1]);
                        write(nounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),"
                                      "' result ',i2,' is',1p,a)"),
			  n, jtype, ioldsd, jr, buf;
                    }
                }
            }
        //
        statement_180:;
        }
    }
    //
    //     Summary
    //
    Alasvm("DGS", nounit, nerrs, ntestt, 0);
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Rdrges
    //
}

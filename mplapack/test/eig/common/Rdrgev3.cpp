/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DDRGEV3.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

void Rdrgev3(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER (&iseed)[4], REAL const thresh, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *b, REAL *s, REAL *t, REAL *q, INTEGER const ldq, REAL *z, REAL *qe, INTEGER const ldqe, REAL *alphar, REAL *alphai, REAL *beta, REAL *alphr1, REAL *alphi1, REAL *beta1, REAL *work, INTEGER const lwork, REAL *result, INTEGER &info) {
    common cmn;
    common_write write(cmn);
    static INTEGER kclass[27] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4};
    static INTEGER kz1[6] = {0, 1, 2, 1, 3, 3};
    static INTEGER kz2[6] = {0, 0, 1, 2, 1, 1};
    static INTEGER kadd[6] = {0, 0, 0, 0, 3, 2};
    static INTEGER katype[27] = {0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0, 0};
    static INTEGER kbtype[27] = {0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0};
    static INTEGER kazero[27] = {1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1, 1};
    static INTEGER kbzero[27] = {1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1, 1};
    static INTEGER kamagn[27] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1, 3};
    static INTEGER kbmagn[27] = {1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1, 3};
    static INTEGER ktrian[27] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    static INTEGER iasign[27] = {0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0, 0};
    static INTEGER ibsign[27] = {0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    REAL safmin = 0.0;
    REAL ulp = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    REAL ulpinv = 0.0;
    REAL rmagn[3 + 1];
    INTEGER ntestt = 0;
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER n1 = 0;
    const INTEGER maxtyp = 27;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ioldsd[4];
    INTEGER ierr = 0;
    INTEGER in = 0;
    INTEGER iadd = 0;
    INTEGER jc = 0;
    INTEGER jr = 0;
    INTEGER i = 0;
    //
    static const char *format_9999 = "(' Rdrgev3: ',a,' returned INFO=',i6,'.',/,3x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',4(i4,','),i5,')')";
    //
    static const char *format_9998 = "(' Rdrgev3: ',a,' Eigenvectors from ',a,' incorrectly normalized.',/,"
                                     "' Bits of error=',0p,g10.3,',',3x,'N=',i4,', JTYPE=',i3,', ISEED=(',4(i4,"
                                     "','),i5,')')";
    //
    static const char *format_9997 = "(/,1x,a3,' -- Real Generalized eigenvalue problem driver')";
    //
    static const char *format_9996 = "(' Matrix types (see Rdrgev3 for details): ')";
    //
    static const char *format_9995 = "(' Special Matrices:',23x,'(J''=transposed Jordan block)',/,"
                                     "'   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ',"
                                     "'6=(diag(J'',I), diag(I,J''))',/,' Diagonal Matrices:  ( ',"
                                     "'D=diag(0,1,2,...) )',/,'   7=(D,I)   9=(large*D, small*I',"
                                     "')  11=(large*I, small*D)  13=(large*D, large*I)',/,"
                                     "'   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ',"
                                     "' 14=(small*D, small*I)',/,'  15=(D, reversed D)')";
    static const char *format_9994 = "(' Matrices Rotated by Random ',a,' Matrices U, V:',/,"
                                     "'  16=Transposed Jordan Blocks             19=geometric ',"
                                     "'alpha, beta=0,1',/,'  17=arithm. alpha&beta             ',"
                                     "'      20=arithmetic alpha, beta=0,1',/,'  18=clustered ',"
                                     "'alpha, beta=0,1            21=random alpha, beta=0,1',/,"
                                     "' Large & Small Matrices:',/,'  22=(large, small)   ',"
                                     "'23=(small,large)    24=(small,small)    25=(large,large)',/,"
                                     "'  26=random O(1) matrices.')";
    //
    static const char *format_9993 = "(/,' Tests performed:    ',/,' 1 = max | ( b A - a B )''*l | / const.,',"
                                     "/,' 2 = | |VR(i)| - 1 | / ulp,',/,"
                                     "' 3 = max | ( b A - a B )*r | / const.',/,' 4 = | |VL(i)| - 1 | / ulp,',"
                                     "/,' 5 = 0 if W same no matter if r or l computed,',/,"
                                     "' 6 = 0 if l same no matter if l computed,',/,"
                                     "' 7 = 0 if r same no matter if r computed,',/,1x)";
    static const char *format_9992 = "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),' result ',i2,"
                                     "' is',0p,f8.2)";
    static const char *format_9991 = "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),' result ',i2,"
                                     "' is',1p,d10.3)";
    //
    // Check for errors
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
    // Compute workspace
    // (Note: Comments in the code beginning "Workspace:" describe the
    // minimal amount of workspace needed at that point in the code,
    // as well as the preferred amount for good performance.
    // NB refers to the optimal block size for the immediately
    // following subroutine, as returned by iMlaenv.
    //
    minwrk = 1;
    if (info == 0 && lwork >= 1) {
        minwrk = max((INTEGER)1, 8 * nmax, nmax * (nmax + 1));
        maxwrk = 7 * nmax + nmax * iMlaenv(1, "Rgeqrf", " ", nmax, 1, nmax, 0);
        maxwrk = max(maxwrk, nmax * (nmax + 1));
        work[1 - 1] = maxwrk;
    }
    //
    if (lwork < minwrk) {
        info = -25;
    }
    //
    if (info != 0) {
        Mxerbla("Rdrgev3", -info);
        return;
    }
    //
    // Quick return if possible
    //
    if (nsizes == 0 || ntypes == 0) {
        return;
    }
    //
    safmin = Rlamch("Safe minimum");
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    safmin = safmin / ulp;
    safmax = one / safmin;
    ulpinv = one / ulp;
    //
    // The values RMAGN(2:3) depend on N, see below.
    //
    rmagn[0] = zero;
    rmagn[1] = one;
    //
    // Loop over sizes, types
    //
    ntestt = 0;
    nerrs = 0;
    nmats = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        n1 = max((INTEGER)1, n);
        rmagn[2] = safmax * ulp / castREAL(n1);
        rmagn[3] = safmin * ulpinv * n1;
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
            // Save ISEED in case of an error.
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            // Generate test matrices A and B
            //
            // Description of control parameters:
            //
            // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
            // =3 means random, =4 means random generalized
            // upper Hessenberg.
            // KATYPE: the "type" to be passed to Rlatm4 for computing A.
            // KAZERO: the pattern of zeros on the diagonal for A:
            // =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
            // =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
            // =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
            // non-zero entries.)
            // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
            // =2: large, =3: small.
            // IASIGN: 1 if the diagonal elements of A are to be
            // multiplied by a random magnitude 1 number, =2 if
            // randomly chosen diagonal blocks are to be rotated
            // to form 2x2 blocks.
            // KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
            // KTRIAN: =0: don't fill in the upper triangle, =1: do.
            // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
            // RMAGN: used to implement KAMAGN and KBMAGN.
            //
            if (mtypes > maxtyp) {
                goto statement_100;
            }
            ierr = 0;
            if (kclass[jtype - 1] < 3) {
                //
                // Generate A (w/o rotation)
                //
                if (abs(katype[jtype - 1]) == 3) {
                    in = 2 * ((n - 1) / 2) + 1;
                    if (in != n) {
                        Rlaset("Full", n, n, zero, zero, a, lda);
                    }
                } else {
                    in = n;
                }
                Rlatm4(katype[jtype - 1], in, kz1[kazero[jtype - 1] - 1], kz2[kazero[jtype - 1] - 1], iasign[jtype - 1], rmagn[kamagn[jtype - 1]], ulp, rmagn[(ktrian[jtype - 1] * kamagn[jtype - 1])], 2, iseed, a, lda);
                iadd = kadd[kazero[jtype - 1] - 1];
                if (iadd > 0 && iadd <= n) {
                    a[(iadd - 1) + (iadd - 1) * lda] = one;
                }
                //
                // Generate B (w/o rotation)
                //
                if (abs(kbtype[jtype - 1]) == 3) {
                    in = 2 * ((n - 1) / 2) + 1;
                    if (in != n) {
                        Rlaset("Full", n, n, zero, zero, b, lda);
                    }
                } else {
                    in = n;
                }
                Rlatm4(kbtype[jtype - 1], in, kz1[kbzero[jtype - 1] - 1], kz2[kbzero[jtype - 1] - 1], ibsign[jtype - 1], rmagn[kbmagn[jtype - 1]], one, rmagn[(ktrian[jtype - 1] * kbmagn[jtype - 1])], 2, iseed, b, lda);
                iadd = kadd[kbzero[jtype - 1] - 1];
                if (iadd != 0 && iadd <= n) {
                    b[(iadd - 1) + (iadd - 1) * lda] = one;
                }
                //
                if (kclass[jtype - 1] == 2 && n > 0) {
                    //
                    // Include rotations
                    //
                    // Generate Q, Z as Householder transformations times
                    // a diagonal matrix.
                    //
                    for (jc = 1; jc <= n - 1; jc = jc + 1) {
                        for (jr = jc; jr <= n; jr = jr + 1) {
                            q[(jr - 1) + (jc - 1) * ldq] = Rlarnd(3, iseed);
                            z[(jr - 1) + (jc - 1) * ldq] = Rlarnd(3, iseed);
                        }
                        Rlarfg(n + 1 - jc, q[(jc - 1) + (jc - 1) * ldq], &q[((jc + 1) - 1) + (jc - 1) * ldq], 1, work[jc - 1]);
                        work[(2 * n + jc) - 1] = sign(one, q[(jc - 1) + (jc - 1) * ldq]);
                        q[(jc - 1) + (jc - 1) * ldq] = one;
                        Rlarfg(n + 1 - jc, z[(jc - 1) + (jc - 1) * ldq], &z[((jc + 1) - 1) + (jc - 1) * ldq], 1, work[(n + jc) - 1]);
                        work[(3 * n + jc) - 1] = sign(one, z[(jc - 1) + (jc - 1) * ldq]);
                        z[(jc - 1) + (jc - 1) * ldq] = one;
                    }
                    q[(n - 1) + (n - 1) * ldq] = one;
                    work[n - 1] = zero;
                    work[(3 * n) - 1] = sign(one, Rlarnd(2, iseed));
                    z[(n - 1) + (n - 1) * ldq] = one;
                    work[(2 * n) - 1] = zero;
                    work[(4 * n) - 1] = sign(one, Rlarnd(2, iseed));
                    //
                    // Apply the diagonal matrices
                    //
                    for (jc = 1; jc <= n; jc = jc + 1) {
                        for (jr = 1; jr <= n; jr = jr + 1) {
                            a[(jr - 1) + (jc - 1) * lda] = work[(2 * n + jr) - 1] * work[(3 * n + jc) - 1] * a[(jr - 1) + (jc - 1) * lda];
                            b[(jr - 1) + (jc - 1) * lda] = work[(2 * n + jr) - 1] * work[(3 * n + jc) - 1] * b[(jr - 1) + (jc - 1) * lda];
                        }
                    }
                    Rorm2r("L", "N", n, n, n - 1, q, ldq, work, a, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                    Rorm2r("R", "T", n, n, n - 1, z, ldq, &work[(n + 1) - 1], a, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                    Rorm2r("L", "N", n, n, n - 1, q, ldq, work, b, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                    Rorm2r("R", "T", n, n, n - 1, z, ldq, &work[(n + 1) - 1], b, lda, &work[(2 * n + 1) - 1], ierr);
                    if (ierr != 0) {
                        goto statement_90;
                    }
                }
            } else if (kclass[jtype - 1] == 3) {
                //
                // Random matrices
                //
                for (jc = 1; jc <= n; jc = jc + 1) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        a[(jr - 1) + (jc - 1) * lda] = rmagn[kamagn[jtype - 1]] * Rlarnd(2, iseed);
                        b[(jr - 1) + (jc - 1) * lda] = rmagn[kbmagn[jtype - 1]] * Rlarnd(2, iseed);
                    }
                }
            } else {
                //
                // Random upper Hessenberg pencil with singular B
                //
                for (jc = 1; jc <= n; jc = jc + 1) {
                    for (jr = 1; jr <= min(jc + 1, n); jr = jr + 1) {
                        a[(jr - 1) + (jc - 1) * lda] = rmagn[kamagn[jtype - 1]] * Rlarnd(2, iseed);
                    }
                    for (jr = jc + 2; jr <= n; jr = jr + 1) {
                        a[(jr - 1) + (jc - 1) * lda] = zero;
                    }
                }
                for (jc = 1; jc <= n; jc = jc + 1) {
                    for (jr = 1; jr <= jc; jr = jr + 1) {
                        b[(jr - 1) + (jc - 1) * lda] = rmagn[kamagn[jtype - 1]] * Rlarnd(2, iseed);
                    }
                    for (jr = jc + 1; jr <= n; jr = jr + 1) {
                        b[(jr - 1) + (jc - 1) * lda] = zero;
                    }
                }
                for (jc = 1; jc <= n; jc = jc + 4) {
                    b[(jc - 1) + (jc - 1) * lda] = zero;
                }
                //
            }
        //
        statement_90:
            //
            if (ierr != 0) {
                write(nounit, format_9999), "Generator", ierr, n, jtype, ioldsd;
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
            // Call Mxlaenv to set the parameters used in Rlaqz0
            //
            Mxlaenv(12, 10);
            Mxlaenv(13, 12);
            Mxlaenv(14, 13);
            Mxlaenv(15, 2);
            Mxlaenv(17, 10);
            //
            // Call Rggev3 to compute eigenvalues and eigenvectors.
            //
            Rlacpy(" ", n, n, a, lda, s, lda);
            Rlacpy(" ", n, n, b, lda, t, lda);
            Rggev3("V", "V", n, s, lda, t, lda, alphar, alphai, beta, q, ldq, z, ldq, work, lwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Rggev31", ierr, n, jtype, ioldsd;
                info = abs(ierr);
                goto statement_190;
            }
            //
            // Do the tests (1) and (2)
            //
            Rget52(true, n, a, lda, b, lda, q, ldq, alphar, alphai, beta, work, &result[1 - 1]);
            if (result[2 - 1] > thresh) {
                write(nounit, format_9998), "Left", "Rggev31", result[2 - 1], n, jtype, ioldsd;
            }
            //
            // Do the tests (3) and (4)
            //
            Rget52(false, n, a, lda, b, lda, z, ldq, alphar, alphai, beta, work, &result[3 - 1]);
            if (result[4 - 1] > thresh) {
                write(nounit, format_9998), "Right", "Rggev31", result[4 - 1], n, jtype, ioldsd;
            }
            //
            // Do the test (5)
            //
            Rlacpy(" ", n, n, a, lda, s, lda);
            Rlacpy(" ", n, n, b, lda, t, lda);
            Rggev3("N", "N", n, s, lda, t, lda, alphr1, alphi1, beta1, q, ldq, z, ldq, work, lwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Rggev32", ierr, n, jtype, ioldsd;
                info = abs(ierr);
                goto statement_190;
            }
            //
            for (j = 1; j <= n; j = j + 1) {
#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
                REAL scale = max(safmin, abs(alphar[j - 1]) + abs(alphai[j - 1]) + abs(beta[j - 1]), abs(alphr1[j - 1]) + abs(alphi1[j - 1]) + abs(beta1[j - 1]));
                REAL diff = abs(alphar[j - 1] - alphr1[j - 1]) + abs(alphai[j - 1] - alphi1[j - 1]) + abs(beta[j - 1] - beta1[j - 1]);
                if (diff > 100.0 * ulp * scale) {
#else
                if (alphar[j - 1] != alphr1[j - 1] || alphai[j - 1] != alphi1[j - 1] || beta[j - 1] != beta1[j - 1]) {
#endif
                    result[5 - 1] = ulpinv;
                }
            }
            //
            // Do the test (6): Compute eigenvalues and left eigenvectors,
            // and test them
            //
            Rlacpy(" ", n, n, a, lda, s, lda);
            Rlacpy(" ", n, n, b, lda, t, lda);
            Rggev3("V", "N", n, s, lda, t, lda, alphr1, alphi1, beta1, qe, ldqe, z, ldq, work, lwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Rggev33", ierr, n, jtype, ioldsd;
                info = abs(ierr);
                goto statement_190;
            }
            //
            for (j = 1; j <= n; j = j + 1) {
#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
                REAL scale = max(safmin, abs(alphar[j - 1]) + abs(alphai[j - 1]) + abs(beta[j - 1]), abs(alphr1[j - 1]) + abs(alphi1[j - 1]) + abs(beta1[j - 1]));
                REAL diff = abs(alphar[j - 1] - alphr1[j - 1]) + abs(alphai[j - 1] - alphi1[j - 1]) + abs(beta[j - 1] - beta1[j - 1]);
                if (diff > 100.0 * ulp * scale) {
#else
                if (alphar[j - 1] != alphr1[j - 1] || alphai[j - 1] != alphi1[j - 1] || beta[j - 1] != beta1[j - 1]) {
#endif
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
            // DO the test (7): Compute eigenvalues and right eigenvectors,
            // and test them
            //
            Rlacpy(" ", n, n, a, lda, s, lda);
            Rlacpy(" ", n, n, b, lda, t, lda);
            Rggev3("N", "V", n, s, lda, t, lda, alphr1, alphi1, beta1, q, ldq, qe, ldqe, work, lwork, ierr);
            if (ierr != 0 && ierr != n + 1) {
                result[1 - 1] = ulpinv;
                write(nounit, format_9999), "Rggev34", ierr, n, jtype, ioldsd;
                info = abs(ierr);
                goto statement_190;
            }
            //
            for (j = 1; j <= n; j = j + 1) {
#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
                REAL scale = max(safmin, abs(alphar[j - 1]) + abs(alphai[j - 1]) + abs(beta[j - 1]), abs(alphr1[j - 1]) + abs(alphi1[j - 1]) + abs(beta1[j - 1]));
                REAL diff = abs(alphar[j - 1] - alphr1[j - 1]) + abs(alphai[j - 1] - alphi1[j - 1]) + abs(beta[j - 1] - beta1[j - 1]);
                if (diff > 100.0 * ulp * scale) {
#else
                if (alphar[j - 1] != alphr1[j - 1] || alphai[j - 1] != alphi1[j - 1] || beta[j - 1] != beta1[j - 1]) {
#endif
                    result[7 - 1] = ulpinv;
                }
            }
            //
            for (j = 1; j <= n; j = j + 1) {
                for (jc = 1; jc <= n; jc = jc + 1) {
                    if (z[(j - 1) + (jc - 1) * ldq] != qe[(j - 1) + (jc - 1) * ldqe]) {
                        result[7 - 1] = ulpinv;
                    }
                }
            }
        //
        // End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_190:
            //
            ntestt += 7;
            //
            // Print out tests which fail.
            //
            for (jr = 1; jr <= 7; jr = jr + 1) {
                if (result[jr - 1] >= thresh) {
                    //
                    // If this is the first test to fail,
                    // print a header to the data file.
                    //
                    if (nerrs == 0) {
                        write(nounit, format_9997), "DGV";
                        //
                        // Matrix types
                        //
                        write(nounit, format_9996);
                        write(nounit, format_9995);
                        write(nounit, format_9994), "Orthogonal";
                        //
                        // Tests performed
                        //
                        write(nounit, format_9993);
                        //
                    }
                    nerrs++;
                    if (result[jr - 1] < 10000.0) {
                        write(nounit, format_9992), n, jtype, ioldsd, jr, result[jr - 1];
                    } else {
                        write(nounit, format_9991), n, jtype, ioldsd, jr, result[jr - 1];
                    }
                }
            }
        //
        statement_210:;
        }
    }
    //
    // Summary
    //
    Alasvm("DGV", nounit, nerrs, ntestt, 0);
    //
    work[1 - 1] = maxwrk;
    //
    // End of Rdrgev3
    //
}

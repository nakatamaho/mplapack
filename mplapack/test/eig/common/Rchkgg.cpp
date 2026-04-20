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

// Derived from LAPACK routine DCHKGG.
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

void Rchkgg(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER (&iseed)[4], REAL const thresh, bool const tstdif, REAL const thrshn, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *b, REAL *h, REAL *t, REAL *s1, REAL *s2, REAL *p1, REAL *p2, REAL *u, INTEGER const ldu, REAL *v, REAL *q, REAL *z, REAL *alphr1, REAL *alphi1, REAL *beta1, REAL *alphr3, REAL *alphi3, REAL *beta3, REAL *evectl, REAL *evectr, REAL *work, INTEGER const lwork, bool *llwork, REAL *result, INTEGER &info) {
    common cmn;
    common_write write(cmn);
    static INTEGER kclass[26] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3};
    static INTEGER kz1[6] = {0, 1, 2, 1, 3, 3};
    static INTEGER kz2[6] = {0, 0, 1, 2, 1, 1};
    static INTEGER kadd[6] = {0, 0, 0, 0, 3, 2};
    static INTEGER katype[26] = {0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0};
    static INTEGER kbtype[26] = {0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0};
    static INTEGER kazero[26] = {1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1};
    static INTEGER kbzero[26] = {1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1};
    static INTEGER kamagn[26] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1};
    static INTEGER kbmagn[26] = {1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1};
    static INTEGER ktrian[26] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    static INTEGER iasign[26] = {0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2, 0};
    static INTEGER ibsign[26] = {0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
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
    REAL rmagn[3 + 1];
    INTEGER ntestt = 0;
    INTEGER nerrs = 0;
    INTEGER nmats = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER n1 = 0;
    const INTEGER maxtyp = 26;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ntest = 0;
    INTEGER ioldsd[4];
    INTEGER iinfo = 0;
    INTEGER in = 0;
    INTEGER iadd = 0;
    INTEGER jc = 0;
    INTEGER jr = 0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER i1 = 0;
    REAL dumma[4];
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    //
    static const char *format_9999 = "(' Rchkgg: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',3(i5,','),i5,')')";
    //
    static const char *format_9998 = "(' Rchkgg: ',a,' Eigenvectors from ',a,' incorrectly ','normalized.',/,"
                                     "' Bits of error=',0p,g10.3,',',9x,'N=',i6,', JTYPE=',i6,', ISEED=(',3(i5,"
                                     "','),i5,')')";
    //
    static const char *format_9997 = "(/,1x,a3,' -- Real Generalized eigenvalue problem')";
    //
    static const char *format_9996 = "(' Matrix types (see Rchkgg for details): ')";
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
    static const char *format_9993 = "(/,' Tests performed:   (H is Hessenberg, S is Schur, B, ',"
                                     "'T, P are triangular,',/,20x,'U, V, Q, and Z are ',a,', l and r are the',"
                                     "/,20x,'appropriate left and right eigenvectors, resp., a is',/,20x,"
                                     "'alpha, b is beta, and ',a,' means ',a,'.)',/,' 1 = | A - U H V',a,"
                                     "' | / ( |A| n ulp )      2 = | B - U T V',a,' | / ( |B| n ulp )',/,"
                                     "' 3 = | I - UU',a,' | / ( n ulp )             4 = | I - VV',a,"
                                     "' | / ( n ulp )',/,' 5 = | H - Q S Z',a,' | / ( |H| n ulp )',6x,"
                                     "'6 = | T - Q P Z',a,' | / ( |T| n ulp )',/,' 7 = | I - QQ',a,"
                                     "' | / ( n ulp )             8 = | I - ZZ',a,' | / ( n ulp )',/,"
                                     "' 9 = max | ( b S - a P )',a,' l | / const.  10 = max | ( b H - a T )',a,"
                                     "' l | / const.',/,"
                                     "' 11= max | ( b S - a P ) r | / const.   12 = max | ( b H',"
                                     "' - a T ) r | / const.',/,1x)";
    //
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
    // Maximum blocksize and shift -- we assume that blocksize and number
    // of shifts are monotone increasing functions of N.
    //
    lwkopt = max(6 * nmax, 2 * nmax * nmax, (INTEGER)1);
    //
    // Check for errors
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
        Mxerbla("Rchkgg", -info);
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
                goto statement_230;
            }
            nmats++;
            ntest = 0;
            //
            // Save ISEED in case of an error.
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            // Initialize RESULT
            //
            for (j = 1; j <= 15; j = j + 1) {
                result[j - 1] = zero;
            }
            //
            // Compute A and B
            //
            // Description of control parameters:
            //
            // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
            // =3 means random.
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
                goto statement_110;
            }
            iinfo = 0;
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
                    a[(iadd - 1) + (iadd - 1) * lda] = rmagn[kamagn[jtype - 1]];
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
                    b[(iadd - 1) + (iadd - 1) * lda] = rmagn[kbmagn[jtype - 1]];
                }
                //
                if (kclass[jtype - 1] == 2 && n > 0) {
                    //
                    // Include rotations
                    //
                    // Generate U, V as Householder transformations times
                    // a diagonal matrix.
                    //
                    for (jc = 1; jc <= n - 1; jc = jc + 1) {
                        for (jr = jc; jr <= n; jr = jr + 1) {
                            u[(jr - 1) + (jc - 1) * ldu] = Rlarnd(3, iseed);
                            v[(jr - 1) + (jc - 1) * ldu] = Rlarnd(3, iseed);
                        }
                        Rlarfg(n + 1 - jc, u[(jc - 1) + (jc - 1) * ldu], &u[((jc + 1) - 1) + (jc - 1) * ldu], 1, work[jc - 1]);
                        work[(2 * n + jc) - 1] = sign(one, u[(jc - 1) + (jc - 1) * ldu]);
                        u[(jc - 1) + (jc - 1) * ldu] = one;
                        Rlarfg(n + 1 - jc, v[(jc - 1) + (jc - 1) * ldu], &v[((jc + 1) - 1) + (jc - 1) * ldu], 1, work[(n + jc) - 1]);
                        work[(3 * n + jc) - 1] = sign(one, v[(jc - 1) + (jc - 1) * ldu]);
                        v[(jc - 1) + (jc - 1) * ldu] = one;
                    }
                    u[(n - 1) + (n - 1) * ldu] = one;
                    work[n - 1] = zero;
                    work[(3 * n) - 1] = sign(one, Rlarnd(2, iseed));
                    v[(n - 1) + (n - 1) * ldu] = one;
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
                    Rorm2r("L", "N", n, n, n - 1, u, ldu, work, a, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Rorm2r("R", "T", n, n, n - 1, v, ldu, &work[(n + 1) - 1], a, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Rorm2r("L", "N", n, n, n - 1, u, ldu, work, b, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                    Rorm2r("R", "T", n, n, n - 1, v, ldu, &work[(n + 1) - 1], b, lda, &work[(2 * n + 1) - 1], iinfo);
                    if (iinfo != 0) {
                        goto statement_100;
                    }
                }
            } else {
                //
                // Random matrices
                //
                for (jc = 1; jc <= n; jc = jc + 1) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        a[(jr - 1) + (jc - 1) * lda] = rmagn[kamagn[jtype - 1]] * Rlarnd(2, iseed);
                        b[(jr - 1) + (jc - 1) * lda] = rmagn[kbmagn[jtype - 1]] * Rlarnd(2, iseed);
                    }
                }
            }
            //
            anorm = Rlange("1", n, n, a, lda, work);
            bnorm = Rlange("1", n, n, b, lda, work);
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
            // Call Rgeqr2, Rorm2r, and Rgghrd to compute H, T, U, and V
            //
            Rlacpy(" ", n, n, a, lda, h, lda);
            Rlacpy(" ", n, n, b, lda, t, lda);
            ntest = 1;
            result[1 - 1] = ulpinv;
            //
            Rgeqr2(n, n, t, lda, work, &work[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rgeqr2", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Rorm2r("L", "T", n, n, n, t, lda, work, h, lda, &work[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rorm2r", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Rlaset("Full", n, n, zero, one, u, ldu);
            Rorm2r("R", "N", n, n, n, t, lda, work, u, ldu, &work[(n + 1) - 1], iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rorm2r", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Rgghrd("V", "I", n, 1, n, h, lda, t, lda, u, ldu, v, ldu, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rgghrd", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            ntest = 4;
            //
            // Do tests 1--4
            //
            Rget51(1, n, a, lda, h, lda, u, ldu, v, ldu, work, result[1 - 1]);
            Rget51(1, n, b, lda, t, lda, u, ldu, v, ldu, work, result[2 - 1]);
            Rget51(3, n, b, lda, t, lda, u, ldu, u, ldu, work, result[3 - 1]);
            Rget51(3, n, b, lda, t, lda, v, ldu, v, ldu, work, result[4 - 1]);
            //
            // Call Rhgeqz to compute S1, P1, S2, P2, Q, and Z, do tests.
            //
            // Compute T1 and UZ
            //
            // Eigenvalues only
            //
            Rlacpy(" ", n, n, h, lda, s2, lda);
            Rlacpy(" ", n, n, t, lda, p2, lda);
            ntest = 5;
            result[5 - 1] = ulpinv;
            //
            Rhgeqz("E", "N", "N", n, 1, n, s2, lda, p2, lda, alphr3, alphi3, beta3, q, ldu, z, ldu, work, lwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rhgeqz(E)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            // Eigenvalues and Full Schur Form
            //
            Rlacpy(" ", n, n, h, lda, s2, lda);
            Rlacpy(" ", n, n, t, lda, p2, lda);
            //
            Rhgeqz("S", "N", "N", n, 1, n, s2, lda, p2, lda, alphr1, alphi1, beta1, q, ldu, z, ldu, work, lwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rhgeqz(S)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            // Eigenvalues, Schur Form, and Schur Vectors
            //
            Rlacpy(" ", n, n, h, lda, s1, lda);
            Rlacpy(" ", n, n, t, lda, p1, lda);
            //
            Rhgeqz("S", "I", "I", n, 1, n, s1, lda, p1, lda, alphr1, alphi1, beta1, q, ldu, z, ldu, work, lwork, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rhgeqz(V)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            ntest = 8;
            //
            // Do Tests 5--8
            //
            Rget51(1, n, h, lda, s1, lda, q, ldu, z, ldu, work, result[5 - 1]);
            Rget51(1, n, t, lda, p1, lda, q, ldu, z, ldu, work, result[6 - 1]);
            Rget51(3, n, t, lda, p1, lda, q, ldu, q, ldu, work, result[7 - 1]);
            Rget51(3, n, t, lda, p1, lda, z, ldu, z, ldu, work, result[8 - 1]);
            //
            // Compute the Left and Right Eigenvectors of (S1,P1)
            //
            // 9: Compute the left eigenvector Matrix without
            // back transforming:
            //
            ntest = 9;
            result[9 - 1] = ulpinv;
            //
            // To test "SELECT" option, compute half of the eigenvectors
            // in one call, and half in another
            //
            i1 = n / 2;
            for (j = 1; j <= i1; j = j + 1) {
                llwork[j - 1] = true;
            }
            for (j = i1 + 1; j <= n; j = j + 1) {
                llwork[j - 1] = false;
            }
            //
            Rtgevc("L", "S", llwork, n, s1, lda, p1, lda, evectl, ldu, dumma, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtgevc(L,S1)", iinfo, n, jtype, ioldsd;
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
            Rtgevc("L", "S", llwork, n, s1, lda, p1, lda, &evectl[((i1 + 1) - 1) * ldu], ldu, dumma, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtgevc(L,S2)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Rget52(true, n, s1, lda, p1, lda, evectl, ldu, alphr1, alphi1, beta1, work, &dumma[1 - 1]);
            result[9 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thrshn) {
                write(nounit, format_9998), "Left", "Rtgevc(HOWMNY=S)", dumma[2 - 1], n, jtype, ioldsd;
            }
            //
            // 10: Compute the left eigenvector Matrix with
            // back transforming:
            //
            ntest = 10;
            result[10 - 1] = ulpinv;
            Rlacpy("F", n, n, q, ldu, evectl, ldu);
            Rtgevc("L", "B", llwork, n, s1, lda, p1, lda, evectl, ldu, dumma, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtgevc(L,B)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Rget52(true, n, h, lda, t, lda, evectl, ldu, alphr1, alphi1, beta1, work, &dumma[1 - 1]);
            result[10 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thrshn) {
                write(nounit, format_9998), "Left", "Rtgevc(HOWMNY=B)", dumma[2 - 1], n, jtype, ioldsd;
            }
            //
            // 11: Compute the right eigenvector Matrix without
            // back transforming:
            //
            ntest = 11;
            result[11 - 1] = ulpinv;
            //
            // To test "SELECT" option, compute half of the eigenvectors
            // in one call, and half in another
            //
            i1 = n / 2;
            for (j = 1; j <= i1; j = j + 1) {
                llwork[j - 1] = true;
            }
            for (j = i1 + 1; j <= n; j = j + 1) {
                llwork[j - 1] = false;
            }
            //
            Rtgevc("R", "S", llwork, n, s1, lda, p1, lda, dumma, ldu, evectr, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtgevc(R,S1)", iinfo, n, jtype, ioldsd;
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
            Rtgevc("R", "S", llwork, n, s1, lda, p1, lda, dumma, ldu, &evectr[((i1 + 1) - 1) * ldu], ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtgevc(R,S2)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Rget52(false, n, s1, lda, p1, lda, evectr, ldu, alphr1, alphi1, beta1, work, &dumma[1 - 1]);
            result[11 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thresh) {
                write(nounit, format_9998), "Right", "Rtgevc(HOWMNY=S)", dumma[2 - 1], n, jtype, ioldsd;
            }
            //
            // 12: Compute the right eigenvector Matrix with
            // back transforming:
            //
            ntest = 12;
            result[12 - 1] = ulpinv;
            Rlacpy("F", n, n, z, ldu, evectr, ldu);
            Rtgevc("R", "B", llwork, n, s1, lda, p1, lda, dumma, ldu, evectr, ldu, n, in, work, iinfo);
            if (iinfo != 0) {
                write(nounit, format_9999), "Rtgevc(R,B)", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                goto statement_210;
            }
            //
            Rget52(false, n, h, lda, t, lda, evectr, ldu, alphr1, alphi1, beta1, work, &dumma[1 - 1]);
            result[12 - 1] = dumma[1 - 1];
            if (dumma[2 - 1] > thresh) {
                write(nounit, format_9998), "Right", "Rtgevc(HOWMNY=B)", dumma[2 - 1], n, jtype, ioldsd;
            }
            //
            // Tests 13--15 are done only on request
            //
            if (tstdif) {
                //
                // Do Tests 13--14
                //
                Rget51(2, n, s1, lda, s2, lda, q, ldu, z, ldu, work, result[13 - 1]);
                Rget51(2, n, p1, lda, p2, lda, q, ldu, z, ldu, work, result[14 - 1]);
                //
                // Do Test 15
                //
                temp1 = zero;
                temp2 = zero;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max(temp1, abs(alphr1[j - 1] - alphr3[j - 1]) + abs(alphi1[j - 1] - alphi3[j - 1]));
                    temp2 = max(temp2, abs(beta1[j - 1] - beta3[j - 1]));
                }
                //
                temp1 = temp1 / max(safmin, ulp * max(temp1, anorm));
                temp2 = temp2 / max(safmin, ulp * max(temp2, bnorm));
                result[15 - 1] = max(temp1, temp2);
                ntest = 15;
            } else {
                result[13 - 1] = zero;
                result[14 - 1] = zero;
                result[15 - 1] = zero;
                ntest = 12;
            }
        //
        // End of Loop -- Check for RESULT(j) > THRESH
        //
        statement_210:
            //
            ntestt += ntest;
            //
            // Print out tests which fail.
            //
            for (jr = 1; jr <= ntest; jr = jr + 1) {
                if (result[jr - 1] >= thresh) {
                    //
                    // If this is the first test to fail,
                    // print a header to the data file.
                    //
                    if (nerrs == 0) {
                        write(nounit, format_9997), "DGG";
                        //
                        // Matrix types
                        //
                        write(nounit, format_9996);
                        write(nounit, format_9995);
                        write(nounit, format_9994), "Orthogonal";
                        //
                        // Tests performed
                        //
                        {
                            write_loop wloop(cmn, nounit, format_9993);
                            wloop, "orthogonal", "'", "transpose";
                            for (j = 1; j <= 10; j = j + 1) {
                                wloop, "'";
                            }
                        }
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
        statement_230:;
        }
    }
    //
    // Summary
    //
    Rlasum("DGG", nounit, nerrs, ntestt);
    //
    // End of Rchkgg
    //
}

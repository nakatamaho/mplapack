/*
 * Copyright (c) 2008-2021
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

void Clalsd(const char *uplo, INTEGER const smlsiz, INTEGER const n, INTEGER const nrhs, REAL *d, REAL *e, COMPLEX *b, INTEGER const ldb, REAL const rcond, INTEGER &rank, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER &info) {
    //
    //  -- LAPACK computational routine --
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    //
    if (n < 0) {
        info = -3;
    } else if (nrhs < 1) {
        info = -4;
    } else if ((ldb < 1) || (ldb < n)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Clalsd", -info);
        return;
    }
    //
    REAL eps = Rlamch("Epsilon");
    //
    //     Set up the tolerance.
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL rcnd = 0.0;
    if ((rcond <= zero) || (rcond >= one)) {
        rcnd = eps;
    } else {
        rcnd = rcond;
    }
    //
    rank = 0;
    //
    //     Quick return if possible.
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (n == 0) {
        return;
    } else if (n == 1) {
        if (d[1 - 1] == zero) {
            Claset("A", 1, nrhs, czero, czero, b, ldb);
        } else {
            rank = 1;
            Clascl("G", 0, 0, &d[1 - 1], one, 1, nrhs, b, ldb, info);
            d[1 - 1] = abs(d[1 - 1]);
        }
        return;
    }
    //
    //     Rotate the matrix if it is lower bidiagonal.
    //
    INTEGER i = 0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    REAL r = 0.0;
    INTEGER j = 0;
    if (uplo == "L") {
        for (i = 1; i <= n - 1; i = i + 1) {
            Rlartg(d[i - 1], &e[i - 1], cs, sn, r);
            d[i - 1] = r;
            e[i - 1] = sn * d[(i + 1) - 1];
            d[(i + 1) - 1] = cs * d[(i + 1) - 1];
            if (nrhs == 1) {
                CRrot(1, &b[(i - 1)], 1, &b[((i + 1) - 1)], 1, cs, sn);
            } else {
                rwork[(i * 2 - 1) - 1] = cs;
                rwork[(i * 2) - 1] = sn;
            }
        }
        if (nrhs > 1) {
            for (i = 1; i <= nrhs; i = i + 1) {
                for (j = 1; j <= n - 1; j = j + 1) {
                    cs = rwork[(j * 2 - 1) - 1];
                    sn = rwork[(j * 2) - 1];
                    CRrot(1, &b[(j - 1) + (i - 1) * ldb], 1, &b[((j + 1) - 1) + (i - 1) * ldb], 1, cs, sn);
                }
            }
        }
    }
    //
    //     Scale.
    //
    INTEGER nm1 = n - 1;
    REAL orgnrm = Rlanst("M", n, d, e);
    if (orgnrm == zero) {
        Claset("A", n, nrhs, czero, czero, b, ldb);
        return;
    }
    //
    Rlascl("G", 0, 0, orgnrm, one, n, 1, d, n, info);
    Rlascl("G", 0, 0, orgnrm, one, nm1, 1, e, nm1, info);
    //
    //     If N is smaller than the minimum divide size SMLSIZ, then solve
    //     the problem with another solver.
    //
    INTEGER irwu = 0;
    INTEGER irwvt = 0;
    INTEGER irwwrk = 0;
    INTEGER irwrb = 0;
    INTEGER irwib = 0;
    INTEGER irwb = 0;
    INTEGER jcol = 0;
    INTEGER jrow = 0;
    INTEGER jreal = 0;
    INTEGER jimag = 0;
    REAL tol = 0.0;
    if (n <= smlsiz) {
        irwu = 1;
        irwvt = irwu + n * n;
        irwwrk = irwvt + n * n;
        irwrb = irwwrk;
        irwib = irwrb + n * nrhs;
        irwb = irwib + n * nrhs;
        Rlaset("A", n, n, zero, one, &rwork[irwu - 1], n);
        Rlaset("A", n, n, zero, one, &rwork[irwvt - 1], n);
        Rlasdq("U", 0, n, n, n, 0, d, e, &rwork[irwvt - 1], n, &rwork[irwu - 1], n, &rwork[irwwrk - 1], 1, &rwork[irwwrk - 1], info);
        if (info != 0) {
            return;
        }
        //
        //        In the real version, B is passed to Rlasdq and multiplied
        //        internally by Q**H. Here B is complex and that product is
        //        computed below in two steps (real and imaginary parts).
        //
        j = irwb - 1;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].real();
            }
        }
        Rgemm("T", "N", n, nrhs, n, one, &rwork[irwu - 1], n, &rwork[irwb - 1], n, zero, &rwork[irwrb - 1], n);
        j = irwb - 1;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].imag();
            }
        }
        Rgemm("T", "N", n, nrhs, n, one, &rwork[irwu - 1], n, &rwork[irwb - 1], n, zero, &rwork[irwib - 1], n);
        jreal = irwrb - 1;
        jimag = irwib - 1;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                jreal++;
                jimag++;
                b[(jrow - 1) + (jcol - 1) * ldb] = COMPLEX(rwork[jreal - 1], &rwork[jimag - 1]);
            }
        }
        //
        tol = rcnd * abs(diRamax(n, d, 1));
        for (i = 1; i <= n; i = i + 1) {
            if (d[i - 1] <= tol) {
                Claset("A", 1, nrhs, czero, czero, &b[(i - 1)], ldb);
            } else {
                Clascl("G", 0, 0, &d[i - 1], one, 1, nrhs, &b[(i - 1)], ldb, info);
                rank++;
            }
        }
        //
        //        Since B is complex, the following call to Rgemm is performed
        //        in two steps (real and imaginary parts). That is for V * B
        //        (in the real version of the code V**H is stored in WORK).
        //
        //        CALL Rgemm( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO,
        //    $               WORK( NWORK ), N )
        //
        j = irwb - 1;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].real();
            }
        }
        Rgemm("T", "N", n, nrhs, n, one, &rwork[irwvt - 1], n, &rwork[irwb - 1], n, zero, &rwork[irwrb - 1], n);
        j = irwb - 1;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                j++;
                rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].imag();
            }
        }
        Rgemm("T", "N", n, nrhs, n, one, &rwork[irwvt - 1], n, &rwork[irwb - 1], n, zero, &rwork[irwib - 1], n);
        jreal = irwrb - 1;
        jimag = irwib - 1;
        for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
            for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                jreal++;
                jimag++;
                b[(jrow - 1) + (jcol - 1) * ldb] = COMPLEX(rwork[jreal - 1], &rwork[jimag - 1]);
            }
        }
        //
        //        Unscale.
        //
        Rlascl("G", 0, 0, one, orgnrm, n, 1, d, n, info);
        Rlasrt("D", n, d, info);
        Clascl("G", 0, 0, orgnrm, one, n, nrhs, b, ldb, info);
        //
        return;
    }
    //
    //     Book-keeping and setting up some constants.
    //
    const REAL two = 2.0;
    INTEGER nlvl = int(log(n.real() / (smlsiz + 1).real()) / log(two)) + 1;
    //
    INTEGER smlszp = smlsiz + 1;
    //
    INTEGER u = 1;
    INTEGER vt = 1 + smlsiz * n;
    INTEGER difl = vt + smlszp * n;
    INTEGER difr = difl + nlvl * n;
    INTEGER z = difr + nlvl * n * 2;
    INTEGER c = z + nlvl * n;
    INTEGER s = c + n;
    INTEGER poles = s + n;
    INTEGER givnum = poles + 2 * nlvl * n;
    INTEGER nrwork = givnum + 2 * nlvl * n;
    INTEGER bx = 1;
    //
    irwrb = nrwork;
    irwib = irwrb + smlsiz * nrhs;
    irwb = irwib + smlsiz * nrhs;
    //
    INTEGER sizei = 1 + n;
    INTEGER k = sizei + n;
    INTEGER givptr = k + n;
    INTEGER perm = givptr + n;
    INTEGER givcol = perm + nlvl * n;
    INTEGER iwk = givcol + nlvl * n * 2;
    //
    INTEGER st = 1;
    INTEGER sqre = 0;
    INTEGER icmpq1 = 1;
    INTEGER icmpq2 = 0;
    INTEGER nsub = 0;
    //
    for (i = 1; i <= n; i = i + 1) {
        if (abs(d[i - 1]) < eps) {
            d[i - 1] = sign(eps, &d[i - 1]);
        }
    }
    //
    INTEGER nsize = 0;
    INTEGER st1 = 0;
    INTEGER bxst = 0;
    for (i = 1; i <= nm1; i = i + 1) {
        if ((abs(e[i - 1]) < eps) || (i == nm1)) {
            nsub++;
            iwork[nsub - 1] = st;
            //
            //           Subproblem found. First determine its size and then
            //           apply divide and conquer on it.
            //
            if (i < nm1) {
                //
                //              A subproblem with E(I) small for I < NM1.
                //
                nsize = i - st + 1;
                iwork[(sizei + nsub - 1) - 1] = nsize;
            } else if (abs(e[i - 1]) >= eps) {
                //
                //              A subproblem with E(NM1) not too small but I = NM1.
                //
                nsize = n - st + 1;
                iwork[(sizei + nsub - 1) - 1] = nsize;
            } else {
                //
                //              A subproblem with E(NM1) small. This implies an
                //              1-by-1 subproblem at D(N), which is not solved
                //              explicitly.
                //
                nsize = i - st + 1;
                iwork[(sizei + nsub - 1) - 1] = nsize;
                nsub++;
                iwork[nsub - 1] = n;
                iwork[(sizei + nsub - 1) - 1] = 1;
                Ccopy(nrhs, &b[(n - 1)], ldb, &work[(bx + nm1) - 1], n);
            }
            st1 = st - 1;
            if (nsize == 1) {
                //
                //              This is a 1-by-1 subproblem and is not solved
                //              explicitly.
                //
                Ccopy(nrhs, &b[(st - 1)], ldb, &work[(bx + st1) - 1], n);
            } else if (nsize <= smlsiz) {
                //
                //              This is a small subproblem and is solved by Rlasdq.
                //
                Rlaset("A", nsize, nsize, zero, one, &rwork[(vt + st1) - 1], n);
                Rlaset("A", nsize, nsize, zero, one, &rwork[(u + st1) - 1], n);
                Rlasdq("U", 0, nsize, nsize, nsize, 0, &d[st - 1], &e[st - 1], &rwork[(vt + st1) - 1], n, &rwork[(u + st1) - 1], n, &rwork[nrwork - 1], 1, &rwork[nrwork - 1], info);
                if (info != 0) {
                    return;
                }
                //
                //              In the real version, B is passed to Rlasdq and multiplied
                //              internally by Q**H. Here B is complex and that product is
                //              computed below in two steps (real and imaginary parts).
                //
                j = irwb - 1;
                for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
                    for (jrow = st; jrow <= st + nsize - 1; jrow = jrow + 1) {
                        j++;
                        rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].real();
                    }
                }
                Rgemm("T", "N", nsize, nrhs, nsize, one, &rwork[(u + st1) - 1], n, &rwork[irwb - 1], nsize, zero, &rwork[irwrb - 1], nsize);
                j = irwb - 1;
                for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
                    for (jrow = st; jrow <= st + nsize - 1; jrow = jrow + 1) {
                        j++;
                        rwork[j - 1] = b[(jrow - 1) + (jcol - 1) * ldb].imag();
                    }
                }
                Rgemm("T", "N", nsize, nrhs, nsize, one, &rwork[(u + st1) - 1], n, &rwork[irwb - 1], nsize, zero, &rwork[irwib - 1], nsize);
                jreal = irwrb - 1;
                jimag = irwib - 1;
                for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
                    for (jrow = st; jrow <= st + nsize - 1; jrow = jrow + 1) {
                        jreal++;
                        jimag++;
                        b[(jrow - 1) + (jcol - 1) * ldb] = COMPLEX(rwork[jreal - 1], &rwork[jimag - 1]);
                    }
                }
                //
                Clacpy("A", nsize, nrhs, &b[(st - 1)], ldb, &work[(bx + st1) - 1], n);
            } else {
                //
                //              A large problem. Solve it using divide and conquer.
                //
                Rlasda(icmpq1, smlsiz, nsize, sqre, &d[st - 1], &e[st - 1], &rwork[(u + st1) - 1], n, &rwork[(vt + st1) - 1], &iwork[(k + st1) - 1], &rwork[(difl + st1) - 1], &rwork[(difr + st1) - 1], &rwork[(z + st1) - 1], &rwork[(poles + st1) - 1], &iwork[(givptr + st1) - 1], &iwork[(givcol + st1) - 1], n, &iwork[(perm + st1) - 1], &rwork[(givnum + st1) - 1], &rwork[(c + st1) - 1], &rwork[(s + st1) - 1], &rwork[nrwork - 1], &iwork[iwk - 1], info);
                if (info != 0) {
                    return;
                }
                bxst = bx + st1;
                Clalsa(icmpq2, smlsiz, nsize, nrhs, &b[(st - 1)], ldb, &work[bxst - 1], n, &rwork[(u + st1) - 1], n, &rwork[(vt + st1) - 1], &iwork[(k + st1) - 1], &rwork[(difl + st1) - 1], &rwork[(difr + st1) - 1], &rwork[(z + st1) - 1], &rwork[(poles + st1) - 1], &iwork[(givptr + st1) - 1], &iwork[(givcol + st1) - 1], n, &iwork[(perm + st1) - 1], &rwork[(givnum + st1) - 1], &rwork[(c + st1) - 1], &rwork[(s + st1) - 1], &rwork[nrwork - 1], &iwork[iwk - 1], info);
                if (info != 0) {
                    return;
                }
            }
            st = i + 1;
        }
    }
    //
    //     Apply the singular values and treat the tiny ones as zero.
    //
    tol = rcnd * abs(diRamax(n, d, 1));
    //
    for (i = 1; i <= n; i = i + 1) {
        //
        //        Some of the elements in D can be negative because 1-by-1
        //        subproblems were not solved explicitly.
        //
        if (abs(d[i - 1]) <= tol) {
            Claset("A", 1, nrhs, czero, czero, &work[(bx + i - 1) - 1], n);
        } else {
            rank++;
            Clascl("G", 0, 0, &d[i - 1], one, 1, nrhs, &work[(bx + i - 1) - 1], n, info);
        }
        d[i - 1] = abs(d[i - 1]);
    }
    //
    //     Now apply back the right singular vectors.
    //
    icmpq2 = 1;
    for (i = 1; i <= nsub; i = i + 1) {
        st = iwork[i - 1];
        st1 = st - 1;
        nsize = iwork[(sizei + i - 1) - 1];
        bxst = bx + st1;
        if (nsize == 1) {
            Ccopy(nrhs, &work[bxst - 1], n, &b[(st - 1)], ldb);
        } else if (nsize <= smlsiz) {
            //
            //           Since B and BX are complex, the following call to Rgemm
            //           is performed in two steps (real and imaginary parts).
            //
            //           CALL Rgemm( 'T', 'N', NSIZE, NRHS, NSIZE, ONE,
            //    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO,
            //    $                  B( ST, 1 ), LDB )
            //
            j = bxst - n - 1;
            jreal = irwb - 1;
            for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
                j += n;
                for (jrow = 1; jrow <= nsize; jrow = jrow + 1) {
                    jreal++;
                    rwork[jreal - 1] = work[(j + jrow) - 1].real();
                }
            }
            Rgemm("T", "N", nsize, nrhs, nsize, one, &rwork[(vt + st1) - 1], n, &rwork[irwb - 1], nsize, zero, &rwork[irwrb - 1], nsize);
            j = bxst - n - 1;
            jimag = irwb - 1;
            for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
                j += n;
                for (jrow = 1; jrow <= nsize; jrow = jrow + 1) {
                    jimag++;
                    rwork[jimag - 1] = work[(j + jrow) - 1].imag();
                }
            }
            Rgemm("T", "N", nsize, nrhs, nsize, one, &rwork[(vt + st1) - 1], n, &rwork[irwb - 1], nsize, zero, &rwork[irwib - 1], nsize);
            jreal = irwrb - 1;
            jimag = irwib - 1;
            for (jcol = 1; jcol <= nrhs; jcol = jcol + 1) {
                for (jrow = st; jrow <= st + nsize - 1; jrow = jrow + 1) {
                    jreal++;
                    jimag++;
                    b[(jrow - 1) + (jcol - 1) * ldb] = COMPLEX(rwork[jreal - 1], &rwork[jimag - 1]);
                }
            }
        } else {
            Clalsa(icmpq2, smlsiz, nsize, nrhs, &work[bxst - 1], n, &b[(st - 1)], ldb, &rwork[(u + st1) - 1], n, &rwork[(vt + st1) - 1], &iwork[(k + st1) - 1], &rwork[(difl + st1) - 1], &rwork[(difr + st1) - 1], &rwork[(z + st1) - 1], &rwork[(poles + st1) - 1], &iwork[(givptr + st1) - 1], &iwork[(givcol + st1) - 1], n, &iwork[(perm + st1) - 1], &rwork[(givnum + st1) - 1], &rwork[(c + st1) - 1], &rwork[(s + st1) - 1], &rwork[nrwork - 1], &iwork[iwk - 1], info);
            if (info != 0) {
                return;
            }
        }
    }
    //
    //     Unscale and sort the singular values.
    //
    Rlascl("G", 0, 0, one, orgnrm, n, 1, d, n, info);
    Rlasrt("D", n, d, info);
    Clascl("G", 0, 0, orgnrm, one, n, nrhs, b, ldb, info);
    //
    //     End of Clalsd
    //
}

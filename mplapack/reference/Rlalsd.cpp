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

void Rlalsd(const char *uplo, INTEGER const smlsiz, INTEGER const n, INTEGER const nrhs, REAL *d, REAL *e, REAL *b, INTEGER const ldb, REAL const rcond, INTEGER &rank, REAL *work, INTEGER *iwork, INTEGER &info) {
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
        Mxerbla("Rlalsd", -info);
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
    if (n == 0) {
        return;
    } else if (n == 1) {
        if (d[1 - 1] == zero) {
            Rlaset("A", 1, nrhs, zero, zero, b, ldb);
        } else {
            rank = 1;
            Rlascl("G", 0, 0, d[1 - 1], one, 1, nrhs, b, ldb, info);
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
    if (Mlsame(uplo, "L")) {
        for (i = 1; i <= n - 1; i = i + 1) {
            Rlartg(d[i - 1], e[i - 1], cs, sn, r);
            d[i - 1] = r;
            e[i - 1] = sn * d[(i + 1) - 1];
            d[(i + 1) - 1] = cs * d[(i + 1) - 1];
            if (nrhs == 1) {
                Rrot(1, &b[(i - 1)], 1, &b[((i + 1) - 1)], 1, cs, sn);
            } else {
                work[(i * 2 - 1) - 1] = cs;
                work[(i * 2) - 1] = sn;
            }
        }
        if (nrhs > 1) {
            for (i = 1; i <= nrhs; i = i + 1) {
                for (j = 1; j <= n - 1; j = j + 1) {
                    cs = work[(j * 2 - 1) - 1];
                    sn = work[(j * 2) - 1];
                    Rrot(1, &b[(j - 1) + (i - 1) * ldb], 1, &b[((j + 1) - 1) + (i - 1) * ldb], 1, cs, sn);
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
        Rlaset("A", n, nrhs, zero, zero, b, ldb);
        return;
    }
    //
    Rlascl("G", 0, 0, orgnrm, one, n, 1, d, n, info);
    Rlascl("G", 0, 0, orgnrm, one, nm1, 1, e, nm1, info);
    //
    //     If N is smaller than the minimum divide size SMLSIZ, then solve
    //     the problem with another solver.
    //
    INTEGER nwork = 0;
    REAL tol = 0.0;
    if (n <= smlsiz) {
        nwork = 1 + n * n;
        Rlaset("A", n, n, zero, one, work, n);
        Rlasdq("U", 0, n, n, 0, nrhs, d, e, work, n, work, n, b, ldb, &work[nwork - 1], info);
        if (info != 0) {
            return;
        }
        tol = rcnd * abs(d[iRamax(n, d, 1) - 1]);
        for (i = 1; i <= n; i = i + 1) {
            if (d[i - 1] <= tol) {
                Rlaset("A", 1, nrhs, zero, zero, &b[(i - 1)], ldb);
            } else {
                Rlascl("G", 0, 0, d[i - 1], one, 1, nrhs, &b[(i - 1)], ldb, info);
                rank++;
            }
        }
        Rgemm("T", "N", n, nrhs, n, one, work, n, b, ldb, zero, &work[nwork - 1], n);
        Rlacpy("A", n, nrhs, &work[nwork - 1], n, b, ldb);
        //
        //        Unscale.
        //
        Rlascl("G", 0, 0, one, orgnrm, n, 1, d, n, info);
        Rlasrt("D", n, d, info);
        Rlascl("G", 0, 0, orgnrm, one, n, nrhs, b, ldb, info);
        //
        return;
    }
    //
    //     Book-keeping and setting up some constants.
    //
    const REAL two = 2.0;
    INTEGER nlvl = castINTEGER(log(castREAL(n) / castREAL(smlsiz + 1)) / log(two)) + 1;
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
    INTEGER bx = givnum + 2 * nlvl * n;
    nwork = bx + n * nrhs;
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
            d[i - 1] = sign(eps, d[i - 1]);
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
                Rcopy(nrhs, &b[(n - 1)], ldb, &work[(bx + nm1) - 1], n);
            }
            st1 = st - 1;
            if (nsize == 1) {
                //
                //              This is a 1-by-1 subproblem and is not solved
                //              explicitly.
                //
                Rcopy(nrhs, &b[(st - 1)], ldb, &work[(bx + st1) - 1], n);
            } else if (nsize <= smlsiz) {
                //
                //              This is a small subproblem and is solved by Rlasdq.
                //
                Rlaset("A", nsize, nsize, zero, one, &work[(vt + st1) - 1], n);
                Rlasdq("U", 0, nsize, nsize, 0, nrhs, &d[st - 1], &e[st - 1], &work[(vt + st1) - 1], n, &work[nwork - 1], n, &b[(st - 1)], ldb, &work[nwork - 1], info);
                if (info != 0) {
                    return;
                }
                Rlacpy("A", nsize, nrhs, &b[(st - 1)], ldb, &work[(bx + st1) - 1], n);
            } else {
                //
                //              A large problem. Solve it using divide and conquer.
                //
                Rlasda(icmpq1, smlsiz, nsize, sqre, &d[st - 1], &e[st - 1], &work[(u + st1) - 1], n, &work[(vt + st1) - 1], &iwork[(k + st1) - 1], &work[(difl + st1) - 1], &work[(difr + st1) - 1], &work[(z + st1) - 1], &work[(poles + st1) - 1], &iwork[(givptr + st1) - 1], &iwork[(givcol + st1) - 1], n, &iwork[(perm + st1) - 1], &work[(givnum + st1) - 1], &work[(c + st1) - 1], &work[(s + st1) - 1], &work[nwork - 1], &iwork[iwk - 1], info);
                if (info != 0) {
                    return;
                }
                bxst = bx + st1;
                Rlalsa(icmpq2, smlsiz, nsize, nrhs, &b[(st - 1)], ldb, &work[bxst - 1], n, &work[(u + st1) - 1], n, &work[(vt + st1) - 1], &iwork[(k + st1) - 1], &work[(difl + st1) - 1], &work[(difr + st1) - 1], &work[(z + st1) - 1], &work[(poles + st1) - 1], &iwork[(givptr + st1) - 1], &iwork[(givcol + st1) - 1], n, &iwork[(perm + st1) - 1], &work[(givnum + st1) - 1], &work[(c + st1) - 1], &work[(s + st1) - 1], &work[nwork - 1], &iwork[iwk - 1], info);
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
    tol = rcnd * abs(d[iRamax(n, d, 1) - 1]);
    //
    for (i = 1; i <= n; i = i + 1) {
        //
        //        Some of the elements in D can be negative because 1-by-1
        //        subproblems were not solved explicitly.
        //
        if (abs(d[i - 1]) <= tol) {
            Rlaset("A", 1, nrhs, zero, zero, &work[(bx + i - 1) - 1], n);
        } else {
            rank++;
            Rlascl("G", 0, 0, d[i - 1], one, 1, nrhs, &work[(bx + i - 1) - 1], n, info);
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
            Rcopy(nrhs, &work[bxst - 1], n, &b[(st - 1)], ldb);
        } else if (nsize <= smlsiz) {
            Rgemm("T", "N", nsize, nrhs, nsize, one, &work[(vt + st1) - 1], n, &work[bxst - 1], n, zero, &b[(st - 1)], ldb);
        } else {
            Rlalsa(icmpq2, smlsiz, nsize, nrhs, &work[bxst - 1], n, &b[(st - 1)], ldb, &work[(u + st1) - 1], n, &work[(vt + st1) - 1], &iwork[(k + st1) - 1], &work[(difl + st1) - 1], &work[(difr + st1) - 1], &work[(z + st1) - 1], &work[(poles + st1) - 1], &iwork[(givptr + st1) - 1], &iwork[(givcol + st1) - 1], n, &iwork[(perm + st1) - 1], &work[(givnum + st1) - 1], &work[(c + st1) - 1], &work[(s + st1) - 1], &work[nwork - 1], &iwork[iwk - 1], info);
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
    Rlascl("G", 0, 0, orgnrm, one, n, nrhs, b, ldb, info);
    //
    //     End of Rlalsd
    //
}

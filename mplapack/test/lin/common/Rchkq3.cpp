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
#include <mplapack_lin.h>

void Rchkq3(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER *nxval, REAL const thresh, REAL *a, REAL *copya, REAL *s, REAL *tau, REAL *work, INTEGER *iwork, INTEGER const nout) {
    common_write write(cmn);
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    if (is_called_first_time) {
        static const INTEGER values[] = {1988, 1989, 1990, 1991};
        data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
    }
    char path[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    REAL eps = 0.0;
    INTEGER im = 0;
    INTEGER m = 0;
    INTEGER lda = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER mnmin = 0;
    INTEGER lwork = 0;
    INTEGER imode = 0;
    const INTEGER ntypes = 6;
    INTEGER mode = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER info = 0;
    INTEGER ilow = 0;
    INTEGER istep = 0;
    INTEGER ihigh = 0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    INTEGER nx = 0;
    INTEGER lw = 0;
    const INTEGER ntests = 3;
    REAL result[ntests];
    INTEGER k = 0;
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize constants and the random number seed.
    //
    path[(1 - 1)] = "Double precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "Q3";
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    eps = Rlamch("Epsilon");
    cmn.infot = 0;
    //
    for (im = 1; im <= nm; im = im + 1) {
        //
        //        Do for each value of M in MVAL.
        //
        m = mval[im - 1];
        lda = max((INTEGER)1, m);
        //
        for (in = 1; in <= nn; in = in + 1) {
            //
            //           Do for each value of N in NVAL.
            //
            n = nval[in - 1];
            mnmin = min(m, n);
            lwork = max({(INTEGER)1, m * max(m, n) + 4 * mnmin + max(m, n), m * n + 2 * mnmin + 4 * n});
            //
            for (imode = 1; imode <= ntypes; imode = imode + 1) {
                if (!dotype[imode - 1]) {
                    goto statement_70;
                }
                //
                //              Do for each type of matrix
                //                 1:  zero matrix
                //                 2:  one small singular value
                //                 3:  geometric distribution of singular values
                //                 4:  first n/2 columns fixed
                //                 5:  last n/2 columns fixed
                //                 6:  every second column fixed
                //
                mode = imode;
                if (imode > 3) {
                    mode = 1;
                }
                //
                //              Generate test matrix of size m by n using
                //              singular value distribution indicated by `mode'.
                //
                for (i = 1; i <= n; i = i + 1) {
                    iwork[i - 1] = 0;
                }
                if (imode == 1) {
                    Rlaset("Full", m, n, zero, zero, copya, lda);
                    for (i = 1; i <= mnmin; i = i + 1) {
                        s[i - 1] = zero;
                    }
                } else {
                    Rlatms(m, n, "Uniform", iseed, "Nonsymm", s, mode, one / eps, one, m, n, "No packing", copya, lda, work, info);
                    if (imode >= 4) {
                        if (imode == 4) {
                            ilow = 1;
                            istep = 1;
                            ihigh = max((INTEGER)1, n / 2);
                        } else if (imode == 5) {
                            ilow = max((INTEGER)1, n / 2);
                            istep = 1;
                            ihigh = n;
                        } else if (imode == 6) {
                            ilow = 1;
                            istep = 2;
                            ihigh = n;
                        }
                        for (i = ilow; i <= ihigh; i = i + istep) {
                            iwork[i - 1] = 1;
                        }
                    }
                    Rlaord("Decreasing", mnmin, s, 1);
                }
                //
                for (inb = 1; inb <= nnb; inb = inb + 1) {
                    //
                    //                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
                    //
                    nb = nbval[inb - 1];
                    nx = nxval[inb - 1];
                    //
                    //                 Get a working copy of COPYA into A and a copy of
                    //                 vector IWORK.
                    //
                    Rlacpy("All", m, n, copya, lda, a, lda);
                    icopy(n, &iwork[1 - 1], 1, &iwork[(n + 1) - 1], 1);
                    //
                    //                 Compute the QR factorization with pivoting of A
                    //
                    lw = max((INTEGER)1, 2 * n + nb * (n + 1));
                    //
                    //                 Compute the QP3 factorization of A
                    //
                    Rgeqp3(m, n, a, lda, &iwork[(n + 1) - 1], tau, work, lw, info);
                    //
                    //                 Compute norm(svd(a) - svd(r))
                    //
                    result[1 - 1] = Rqrt12(m, n, a, lda, s, work, lwork);
                    //
                    //                 Compute norm( A*P - Q*R )
                    //
                    result[2 - 1] = Rqpt01(m, n, mnmin, copya, a, lda, tau, &iwork[(n + 1) - 1], work, lwork);
                    //
                    //                 Compute Q'*Q
                    //
                    result[3 - 1] = Rqrt11(m, mnmin, a, lda, tau, work, lwork);
                    //
                    //                 Print information about the tests that did not pass
                    //                 the threshold.
                    //
                    for (k = 1; k <= ntests; k = k + 1) {
                        if (result[k - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            write(nout, "(1x,a,' M =',i5,', N =',i5,', NB =',i4,', type ',i2,"
                                        "', test ',i2,', ratio =',g12.5)"),
                                "Rgeqp3", m, n, nb, imode, k, result(k);
                            nfail++;
                        }
                    }
                    nrun += ntests;
                    //
                }
            statement_70:;
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rchkq3
    //
}

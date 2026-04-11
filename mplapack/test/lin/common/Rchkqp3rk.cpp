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

// Derived from LAPACK routine DCHKQP3RK.
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
#include <mplapack_lin.h>

void Rchkqp3rk(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, INTEGER const nnb, INTEGER *nbval, INTEGER *nxval, REAL const thresh, REAL *a, REAL *copya, REAL *b, REAL *copyb, REAL *s, REAL *tau, REAL *work, INTEGER *iwork, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    static INTEGER iseedy[4] = {1988, 1989, 1990, 1991};
    //
    static const char *format_9999 = "(1x,a,' M =',i5,', N =',i5,', NRHS =',i5,', KMAX =',i5,', ABSTOL =',"
                                     "g12.5,', RELTOL =',g12.5,', NB =',i4,', NX =',i4,', type ',i2,', test ',"
                                     "i2,', ratio =',g12.5)";
    //
    // Initialize constants and the random number seed.
    //
    fem::str<3> path = "Double precision";
    path(2, 3) = "QK";
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    REAL eps = Rlamch("Epsilon");
    infot = 0;
    //
    INTEGER im = 0;
    INTEGER m = 0;
    INTEGER lda = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER minmn = 0;
    INTEGER lwork = 0;
    INTEGER ins = 0;
    INTEGER nrhs = 0;
    fem::str<1> type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    fem::str<1> dist;
    INTEGER info = 0;
    INTEGER imat = 0;
    const INTEGER ntypes = 19;
    const REAL zero = 0.0;
    INTEGER jb_zero = 0;
    INTEGER nb_zero = 0;
    INTEGER nb_gen = 0;
    INTEGER j_inc = 0;
    INTEGER j_first_nz = 0;
    INTEGER ind_offset_gen = 0;
    INTEGER j = 0;
    INTEGER ind_out = 0;
    INTEGER ind_in = 0;
    INTEGER minmnb_gen = 0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    INTEGER nx = 0;
    INTEGER kmax = 0;
    const INTEGER ntests = 5;
    REAL result[ntests];
    REAL abstol = 0.0;
    REAL reltol = 0.0;
    INTEGER lw = 0;
    INTEGER kfact = 0;
    REAL maxc2nrmk = 0.0;
    REAL relmaxc2nrmk = 0.0;
    REAL dtemp = 0.0;
    const REAL bignum = 1e38;
    INTEGER lwork_mqr = 0;
    const REAL one = 1.0;
    REAL rdummy[1];
    INTEGER t = 0;
    for (im = 1; im <= nm; im = im + 1) {
        //
        // Do for each value of M in MVAL.
        //
        m = mval[im - 1];
        lda = max((INTEGER)1, m);
        //
        for (in = 1; in <= nn; in = in + 1) {
            //
            // Do for each value of N in NVAL.
            //
            n = nval[in - 1];
            minmn = min(m, n);
            lwork = max((INTEGER)1, m * max(m, n) + 4 * minmn + max(m, n), m * n + 2 * minmn + 4 * n);
            //
            for (ins = 1; ins <= nns; ins = ins + 1) {
                nrhs = nsval[ins - 1];
                //
                // Set up parameters with Rlatb4 and generate
                // M-by-NRHS B matrix with Rlatms.
                // IMAT = 14:
                // Random matrix, CNDNUM = 2, NORM = ONE,
                // MODE = 3 (geometric distribution of singular values).
                //
                Rlatb4(path, 14, m, nrhs, type, kl, ku, anorm, mode, cndnum, dist);
                //
                srnamt = "Rlatms";
                Rlatms(m, nrhs, dist, iseed, type, s, mode, cndnum, anorm, kl, ku, "No packing", copyb, lda, work, info);
                //
                // Check error code from Rlatms.
                //
                if (info != 0) {
                    Alaerh(path, "Rlatms", info, 0, " ", m, nrhs, -1, -1, -1, 6, nfail, nerrs, nout);
                    continue;
                }
                //
                for (imat = 1; imat <= ntypes; imat = imat + 1) {
                    //
                    // Do the tests only if DOTYPE( IMAT ) is true.
                    //
                    if (!dotype[imat - 1]) {
                        continue;
                    }
                    //
                    // The type of distribution used to generate the random
                    // eigen-/singular values:
                    // ( 'S' for symmetric distribution ) => UNIFORM( -1, 1 )
                    //
                    // Do for each type of NON-SYMMETRIC matrix:                               CNDNUM                     NORM                                     MODE
                    // 1. Zero matrix
                    // 2. Random, Diagonal, CNDNUM = 2                                        CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 3. Random, Upper triangular, CNDNUM = 2                                CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 4. Random, Lower triangular, CNDNUM = 2                                CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 5. Random, First column is zero, CNDNUM = 2                            CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 6. Random, Last MINMN column is zero, CNDNUM = 2                       CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 7. Random, Last N column is zero, CNDNUM = 2                           CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 8. Random, Middle column in MINMN is zero, CNDNUM = 2                  CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 9. Random, First half of MINMN columns are zero, CNDNUM = 2            CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 10. Random, Last columns are zero starting from MINMN/2+1, CNDNUM = 2   CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 11. Random, Half MINMN columns in the middle are zero starting
                    // from  MINMN/2-(MINMN/2)/2+1, CNDNUM = 2                          CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 12. Random, Odd columns are ZERO, CNDNUM = 2                            CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 13. Random, Even columns are ZERO, CNDNUM = 2                           CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 14. Random, CNDNUM = 2                                                  CNDNUM = 2                      ONE                                      3 ( geometric distribution of singular values )
                    // 15. Random, CNDNUM = sqrt(0.1/EPS)                                      CNDNUM = BADC1 = sqrt(0.1/EPS)  ONE                                      3 ( geometric distribution of singular values )
                    // 16. Random, CNDNUM = 0.1/EPS                                            CNDNUM = BADC2 = 0.1/EPS        ONE                                      3 ( geometric distribution of singular values )
                    // 17. Random, CNDNUM = 0.1/EPS,                                           CNDNUM = BADC2 = 0.1/EPS        ONE                                      2 ( one small singular value, S(N)=1/CNDNUM )
                    // one small singular value S(N)=1/CNDNUM
                    // 18. Random, CNDNUM = 2, scaled near underflow                           CNDNUM = 2                      SMALL = SAFMIN
                    // 19. Random, CNDNUM = 2, scaled near overflow                            CNDNUM = 2                      LARGE = 1.0/( 0.25 * ( SAFMIN / EPS ) )  3 ( geometric distribution of singular values )
                    //
                    if (imat == 1) {
                        //
                        // Matrix 1: Zero matrix
                        //
                        Rlaset("Full", m, n, zero, zero, copya, lda);
                        for (i = 1; i <= minmn; i = i + 1) {
                            s[i - 1] = zero;
                        }
                        //
                    } else if ((imat >= 2 && imat <= 4) || (imat >= 14 && imat <= 19)) {
                        //
                        // Matrices 2-5.
                        //
                        // Set up parameters with Rlatb4 and generate a test
                        // matrix with Rlatms.
                        //
                        Rlatb4(path, imat, m, n, type, kl, ku, anorm, mode, cndnum, dist);
                        //
                        srnamt = "Rlatms";
                        Rlatms(m, n, dist, iseed, type, s, mode, cndnum, anorm, kl, ku, "No packing", copya, lda, work, info);
                        //
                        // Check error code from Rlatms.
                        //
                        if (info != 0) {
                            Alaerh(path, "Rlatms", info, 0, " ", m, n, -1, -1, -1, imat, nfail, nerrs, nout);
                            continue;
                        }
                        //
                        Rlaord("Decreasing", minmn, s, 1);
                        //
                    } else if (minmn >= 2 && imat >= 5 && imat <= 13) {
                        //
                        // Rectangular matrices 5-13 that contain zero columns,
                        // only for matrices MINMN >=2.
                        //
                        // JB_ZERO is the column index of ZERO block.
                        // NB_ZERO is the column block size of ZERO block.
                        // NB_GEN is the column blcok size of the
                        // generated block.
                        // J_INC in the non_zero column index increment
                        // for matrix 12 and 13.
                        // J_FIRS_NZ is the index of the first non-zero
                        // column.
                        //
                        if (imat == 5) {
                            //
                            // First column is zero.
                            //
                            jb_zero = 1;
                            nb_zero = 1;
                            nb_gen = n - nb_zero;
                            //
                        } else if (imat == 6) {
                            //
                            // Last column MINMN is zero.
                            //
                            jb_zero = minmn;
                            nb_zero = 1;
                            nb_gen = n - nb_zero;
                            //
                        } else if (imat == 7) {
                            //
                            // Last column N is zero.
                            //
                            jb_zero = n;
                            nb_zero = 1;
                            nb_gen = n - nb_zero;
                            //
                        } else if (imat == 8) {
                            //
                            // Middle column in MINMN is zero.
                            //
                            jb_zero = minmn / 2 + 1;
                            nb_zero = 1;
                            nb_gen = n - nb_zero;
                            //
                        } else if (imat == 9) {
                            //
                            // First half of MINMN columns is zero.
                            //
                            jb_zero = 1;
                            nb_zero = minmn / 2;
                            nb_gen = n - nb_zero;
                            //
                        } else if (imat == 10) {
                            //
                            // Last columns are zero columns,
                            // starting from (MINMN / 2 + 1) column.
                            //
                            jb_zero = minmn / 2 + 1;
                            nb_zero = n - jb_zero + 1;
                            nb_gen = n - nb_zero;
                            //
                        } else if (imat == 11) {
                            //
                            // Half of the columns in the middle of MINMN
                            // columns is zero, starting from
                            // MINMN/2 - (MINMN/2)/2 + 1 column.
                            //
                            jb_zero = minmn / 2 - (minmn / 2) / 2 + 1;
                            nb_zero = minmn / 2;
                            nb_gen = n - nb_zero;
                            //
                        } else if (imat == 12) {
                            //
                            // Odd-numbered columns are zero,
                            //
                            nb_gen = n / 2;
                            nb_zero = n - nb_gen;
                            j_inc = 2;
                            j_first_nz = 2;
                            //
                        } else if (imat == 13) {
                            //
                            // Even-numbered columns are zero.
                            //
                            nb_zero = n / 2;
                            nb_gen = n - nb_zero;
                            j_inc = 2;
                            j_first_nz = 1;
                            //
                        }
                        //
                        // 1) Set the first NB_ZERO columns in COPYA(1:M,1:N)
                        // to zero.
                        //
                        Rlaset("Full", m, nb_zero, zero, zero, copya, lda);
                        //
                        // 2) Generate an M-by-(N-NB_ZERO) matrix with the
                        // chosen singular value distribution
                        // in COPYA(1:M,NB_ZERO+1:N).
                        //
                        Rlatb4(path, imat, m, nb_gen, type, kl, ku, anorm, mode, cndnum, dist);
                        //
                        srnamt = "Rlatms";
                        //
                        ind_offset_gen = nb_zero * lda;
                        //
                        Rlatms(m, nb_gen, dist, iseed, type, s, mode, cndnum, anorm, kl, ku, "No packing", &copya[(ind_offset_gen + 1) - 1], lda, work, info);
                        //
                        // Check error code from Rlatms.
                        //
                        if (info != 0) {
                            Alaerh(path, "Rlatms", info, 0, " ", m, nb_gen, -1, -1, -1, imat, nfail, nerrs, nout);
                            continue;
                        }
                        //
                        // 3) Swap the gererated colums from the right side
                        // NB_GEN-size block in COPYA into correct column
                        // positions.
                        //
                        if (imat == 6 || imat == 7 || imat == 8 || imat == 10 || imat == 11) {
                            //
                            // Move by swapping the generated columns
                            // from the right NB_GEN-size block from
                            // (NB_ZERO+1:NB_ZERO+JB_ZERO)
                            // into columns (1:JB_ZERO-1).
                            //
                            for (j = 1; j <= jb_zero - 1; j = j + 1) {
                                Rswap(m, &copya[((nb_zero + j - 1) * lda + 1) - 1], 1, &copya[((j - 1) * lda + 1) - 1], 1);
                            }
                            //
                        } else if (imat == 12 || imat == 13) {
                            //
                            // ( IMAT = 12, Odd-numbered ZERO columns. )
                            // Swap the generated columns from the right
                            // NB_GEN-size block into the even zero colums in the
                            // left NB_ZERO-size block.
                            //
                            // ( IMAT = 13, Even-numbered ZERO columns. )
                            // Swap the generated columns from the right
                            // NB_GEN-size block into the odd zero colums in the
                            // left NB_ZERO-size block.
                            //
                            for (j = 1; j <= nb_gen; j = j + 1) {
                                ind_out = (nb_zero + j - 1) * lda + 1;
                                ind_in = (j_inc * (j - 1) + (j_first_nz - 1)) * lda + 1;
                                Rswap(m, &copya[ind_out - 1], 1, &copya[ind_in - 1], 1);
                            }
                            //
                        }
                        //
                        // 5) Order the singular values generated by
                        // DLAMTS in decreasing order and add trailing zeros
                        // that correspond to zero columns.
                        // The total number of singular values is MINMN.
                        //
                        minmnb_gen = min(m, nb_gen);
                        //
                        for (i = minmnb_gen + 1; i <= minmn; i = i + 1) {
                            s[i - 1] = zero;
                        }
                        //
                    } else {
                        //
                        // IF(MINMN.LT.2) skip this size for this matrix type.
                        //
                        continue;
                    }
                    //
                    // Initialize a copy array for a pivot array for Rgeqp3rk.
                    //
                    for (i = 1; i <= n; i = i + 1) {
                        iwork[i - 1] = 0;
                    }
                    //
                    for (inb = 1; inb <= nnb; inb = inb + 1) {
                        //
                        // Do for each pair of values (NB,NX) in NBVAL and NXVAL.
                        //
                        nb = nbval[inb - 1];
                        Mxlaenv(1, nb);
                        nx = nxval[inb - 1];
                        Mxlaenv(3, nx);
                        //
                        // We do MIN(M,N)+1 because we need a test for KMAX > N,
                        // when KMAX is larger than MIN(M,N), KMAX should be
                        // KMAX = MIN(M,N)
                        //
                        for (kmax = 0; kmax <= min(m, n) + 1; kmax = kmax + 1) {
                            //
                            // Get a working copy of COPYA into A( 1:M,1:N ).
                            // Get a working copy of COPYB into A( 1:M, (N+1):NRHS ).
                            // Get a working copy of COPYB into into B( 1:M, 1:NRHS ).
                            // Get a working copy of IWORK(1:N) awith zeroes into
                            // which is going to be used as pivot array IWORK( N+1:2N ).
                            // NOTE: IWORK(2N+1:3N) is going to be used as a WORK array
                            // for the routine.
                            //
                            Rlacpy("All", m, n, copya, lda, a, lda);
                            Rlacpy("All", m, nrhs, copyb, lda, &a[(lda * n + 1) - 1], lda);
                            Rlacpy("All", m, nrhs, copyb, lda, b, lda);
                            icopy(n, &iwork[1 - 1], 1, &iwork[(n + 1) - 1], 1);
                            for (i = 1; i <= ntests; i = i + 1) {
                                result[i - 1] = zero;
                            }
                            //
                            abstol = -1.0;
                            reltol = -1.0;
                            //
                            // Compute the QR factorization with pivoting of A
                            //
                            lw = max((INTEGER)1, max(2 * n + nb * (n + nrhs + 1), 3 * n + nrhs - 1));
                            //
                            // Compute Rgeqp3rk factorization of A.
                            //
                            srnamt = "Rgeqp3rk";
                            Rgeqp3rk(m, n, nrhs, kmax, abstol, reltol, a, lda, kfact, maxc2nrmk, relmaxc2nrmk, &iwork[(n + 1) - 1], tau, work, lw, &iwork[(2 * n + 1) - 1], info);
                            //
                            // Check error code from Rgeqp3rk.
                            //
                            if (info < 0) {
                                Alaerh(path, "Rgeqp3rk", info, 0, " ", m, n, nx, -1, nb, imat, nfail, nerrs, nout);
                            }
                            //
                            // Compute test 1:
                            //
                            // This test in only for the full rank factorization of
                            // the matrix A.
                            //
                            // Array S(1:min(M,N)) contains svd(A) the sigular values
                            // of the original matrix A in decreasing absolute value
                            // order. The test computes svd(R), the vector sigular
                            // values of the upper trapezoid of A(1:M,1:N) that
                            // contains the factor R, in decreasing order. The test
                            // returns the ratio:
                            //
                            // 2-norm(svd(R) - svd(A)) / ( max(M,N) * 2-norm(svd(A)) * EPS )
                            //
                            if (kfact == minmn) {
                                //
                                result[1 - 1] = Rqrt12(m, n, a, lda, s, work, lwork);
                                //
                                nrun++;
                                //
                                // End test 1
                                //
                            }
                            //
                            // Compute test 2:
                            //
                            // The test returns the ratio:
                            //
                            // 1-norm( A*P - Q*R ) / ( max(M,N) * 1-norm(A) * EPS )
                            //
                            result[2 - 1] = Rqpt01(m, n, kfact, copya, a, lda, tau, &iwork[(n + 1) - 1], work, lwork);
                            //
                            // Compute test 3:
                            //
                            // The test returns the ratio:
                            //
                            // 1-norm( Q**T * Q - I ) / ( M * EPS )
                            //
                            result[3 - 1] = Rqrt11(m, kfact, a, lda, tau, work, lwork);
                            //
                            nrun += 2;
                            //
                            // Compute test 4:
                            //
                            // This test is only for the factorizations with the
                            // rank greater than 2.
                            // The elements on the diagonal of R should be non-
                            // increasing.
                            //
                            // The test returns the ratio:
                            //
                            // Returns 1.0D+100 if abs(R(K+1,K+1)) > abs(R(K,K)),
                            // K=1:KFACT-1
                            //
                            if (min(kfact, minmn) >= 2) {
                                //
                                for (j = 1; j <= kfact - 1; j = j + 1) {
                                    //
                                    dtemp = ((abs(a[((j - 1) * lda + j) - 1]) - abs(a[((j)*lda + j + 1) - 1])) / abs(a[1 - 1]));
                                    //
                                    if (dtemp < zero) {
                                        result[4 - 1] = bignum;
                                    }
                                    //
                                }
                                //
                                nrun++;
                                //
                                // End test 4.
                                //
                            }
                            //
                            // Compute test 5:
                            //
                            // This test in only for matrix A with min(M,N) > 0.
                            //
                            // The test returns the ratio:
                            //
                            // 1-norm(Q**T * B - Q**T * B ) /
                            // ( M * EPS )
                            //
                            // (1) Compute B:=Q**T * B in the matrix B.
                            //
                            if (minmn > 0) {
                                //
                                lwork_mqr = max((INTEGER)1, nrhs);
                                Rormqr("Left", "Transpose", m, nrhs, kfact, a, lda, tau, b, lda, work, lwork_mqr, info);
                                //
                                for (i = 1; i <= nrhs; i = i + 1) {
                                    //
                                    // Compare N+J-th column of A and J-column of B.
                                    //
                                    Raxpy(m, -one, &a[((n + i - 1) * lda + 1) - 1], 1, &b[((i - 1) * lda + 1) - 1], 1);
                                }
                                //
                                result[5 - 1] = abs(Rlange("One-norm", m, nrhs, b, lda, rdummy) / (castREAL(m) * Rlamch("Epsilon")));
                                //
                                nrun++;
                                //
                                // End compute test 5.
                                //
                            }
                            //
                            // Print information about the tests that did not
                            // pass the threshold.
                            //
                            for (t = 1; t <= ntests; t = t + 1) {
                                if (result[t - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Alahd(nout, path);
                                    }
                                    write(nout, format_9999), "Rgeqp3rk", m, n, nrhs, kmax, abstol, reltol, nb, nx, imat, t, result[t - 1];
                                    nfail++;
                                }
                            }
                            //
                            // END DO KMAX = 1, MIN(M,N)+1
                            //
                        }
                        //
                        // END DO for INB = 1, NNB
                        //
                    }
                    //
                    // END DO  for IMAT = 1, NTYPES
                    //
                }
                //
                // END DO for INS = 1, NNS
                //
            }
            //
            // END DO for IN = 1, NN
            //
        }
        //
        // END DO for IM = 1, NM
        //
    }
    //
    // Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    // End of Rchkqp3rk
    //
}

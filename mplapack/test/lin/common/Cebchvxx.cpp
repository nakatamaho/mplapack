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

void Cebchvxx(common &cmn, REAL const thresh, const char *path) {
    common_write write(cmn);
    //     .. Scalar Arguments ..
    //
    //     .. Local Scalars ..
    //
    //     .. Local Arrays ..
    //
    //     .. External Functions ..
    //
    //     .. External Subroutines ..
    //
    //     .. Intrinsic Functions ..
    //
    //     .. Statement Functions ..
    //
    //     .. Statement Function Definitions ..
    COMPLEX zdum = 0.0;
    abs1(zdum) = abs(zdum.real()) + abs(zdum.imag());
    //
    //     .. Parameters ..
    //
    //  Create the loop to test out the Hilbert matrices
    //
    char fact = "E";
    char uplo = "U";
    char trans = "N";
    char equed = "N";
    REAL eps = Rlamch("Epsilon");
    INTEGER nfail = 0;
    INTEGER n_aux_tests = 0;
    const INTEGER nmax = 10;
    INTEGER lda = nmax;
    INTEGER ldab = (nmax - 1) + (nmax - 1) + 1;
    INTEGER ldafb = 2 * (nmax - 1) + (nmax - 1) + 1;
    str<2> c2 = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Main loop to test the different Hilbert Matrices.
    //
    bool printed_guide = false;
    //
    INTEGER n = 0;
    const INTEGER nparams = 2;
    arr_1d<nparams, REAL> params(fill0);
    INTEGER kl = 0;
    INTEGER ku = 0;
    INTEGER nrhs = 0;
    REAL m = 0.0;
    arr_2d<nmax, nmax, COMPLEX> a(fill0);
    arr_2d<nmax, nmax, COMPLEX> invhilb(fill0);
    arr_2d<nmax, nmax, COMPLEX> b(fill0);
    arr_1d<nmax * 3 * 5, COMPLEX> work(fill0);
    INTEGER info = 0;
    arr_2d<nmax, nmax, COMPLEX> acopy(fill0);
    INTEGER j = 0;
    INTEGER i = 0;
    arr_2d<(nmax - 1) + (nmax - 1) + 1, nmax, COMPLEX> ab(fill0);
    arr_2d<(nmax - 1) + (nmax - 1) + 1, nmax, COMPLEX> abcopy(fill0);
    arr_2d<nmax, nmax, COMPLEX> af(fill0);
    arr_1d<nmax, int> ipiv(fill0);
    arr_1d<nmax, REAL> s(fill0);
    arr_2d<nmax, nmax, COMPLEX> x(fill0);
    REAL orcond = 0.0;
    REAL rpvgrw = 0.0;
    arr_1d<nmax, REAL> berr(fill0);
    const INTEGER nerrbnd = 3;
    arr_1d<nmax * 3, REAL> errbnd_n(fill0);
    arr_1d<nmax * 3, REAL> errbnd_c(fill0);
    arr_1d<3 * nmax, REAL> rwork(fill0);
    nmax), fill0);
    arr_1d<nmax, REAL> r(fill0);
    arr_1d<nmax, REAL> c(fill0);
    REAL rcond = 0.0;
    arr_2d<nmax, nmax, REAL> diff(fill0);
    REAL rnorm = 0.0;
    REAL rinorm = 0.0;
    REAL sumr = 0.0;
    REAL sumri = 0.0;
    arr_1d<nmax, REAL> rinv(fill0);
    REAL ncond = 0.0;
    REAL condthresh = 0.0;
    REAL errthresh = 0.0;
    INTEGER k = 0;
    REAL normt = 0.0;
    REAL normdif = 0.0;
    REAL cwise_err = 0.0;
    REAL nwise_err = 0.0;
    REAL ccond = 0.0;
    const INTEGER bnd_i = 2;
    REAL nwise_bnd = 0.0;
    REAL cwise_bnd = 0.0;
    const INTEGER cond_i = 3;
    REAL nwise_rcond = 0.0;
    REAL cwise_rcond = 0.0;
    char nguar = char0;
    const INTEGER ntests = 6;
    arr_1d<ntests, REAL> tstrat(fill0);
    char cguar = char0;
    for (n = 1; n <= nmax; n = n + 1) {
        params[1 - 1] = -1;
        params[2 - 1] = -1;
        //
        kl = n - 1;
        ku = n - 1;
        nrhs = n;
        m = max(sqrt(n.real()), 10.0);
        //
        //        Generate the Hilbert matrix, its inverse, and the
        //        right hand side, all scaled by the LCM(1,..,2N-1).
        Clahilb(n, n, a, lda, invhilb, lda, b, lda, work, info, path);
        //
        //        Copy A into ACOPY.
        zlacpy("ALL", n, n, a, nmax, acopy, nmax);
        //
        //        Store A in band format for GB tests
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= kl + ku + 1; i = i + 1) {
                ab[(i - 1) + (j - 1) * ldab] = (0.0, 0.0);
            }
        }
        for (j = 1; j <= n; j = j + 1) {
            for (i = max((INTEGER)1, j - ku); i <= min(n, j + kl); i = i + 1) {
                ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab] = a[(i - 1) + (j - 1) * lda];
            }
        }
        //
        //        Copy AB into ABCOPY.
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= kl + ku + 1; i = i + 1) {
                abcopy[(i - 1) + (j - 1) * ldabcopy] = (0.0, 0.0);
            }
        }
        zlacpy("ALL", kl + ku + 1, n, ab, ldab, abcopy, ldab);
        //
        //        Call Z**SVXX with default PARAMS and N_ERR_BND = 3.
        if (Mlsamen2, c2, "SY") {
            zsysvxx(fact, uplo, n, nrhs, acopy, lda, af, lda, ipiv, equed, s, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else if (Mlsamen2, c2, "PO") {
            zposvxx(fact, uplo, n, nrhs, acopy, lda, af, lda, equed, s, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else if (Mlsamen2, c2, "HE") {
            zhesvxx(fact, uplo, n, nrhs, acopy, lda, af, lda, ipiv, equed, s, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else if (Mlsamen2, c2, "GB") {
            zgbsvxx(fact, trans, n, kl, ku, nrhs, abcopy, ldab, afb, ldafb, ipiv, equed, r, c, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else {
            zgesvxx(fact, trans, n, nrhs, acopy, lda, af, lda, ipiv, equed, r, c, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        }
        //
        n_aux_tests++;
        if (orcond < eps) {
            //        Either factorization failed or the matrix is flagged, and 1 <=
            //        INFO <= N+1. We don't decide based on rcond anymore.
            //            IF (INFO .EQ. 0 .OR. INFO .GT. N+1) THEN
            //               NFAIL = NFAIL + 1
            //               WRITE (*, FMT=8000) N, INFO, ORCOND, RCOND
            //            END IF
        } else {
            //        Either everything succeeded (INFO == 0) or some solution failed
            //        to converge (INFO > N+1).
            if (info > 0 && info <= n + 1) {
                nfail++;
                write(6, "(' Z',a2,'SVXX: N =',i2,', INFO = ',i3,', ORCOND = ',g12.5,"
                         "', real RCOND = ',g12.5)"),
                    c2, n, info, orcond, rcond;
            }
        }
        //
        //        Calculating the difference between Z**SVXX's X and the true X.
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= nrhs; j = j + 1) {
                diff[(i - 1) + (j - 1) * lddiff] = x[(i - 1) + (j - 1) * ldx] - invhilb[(i - 1) + (j - 1) * ldinvhilb];
            }
        }
        //
        //        Calculating the RCOND
        rnorm = 0;
        rinorm = 0;
        if (Mlsamen2, c2, "PO" || Mlsamen2, c2, "SY" || Mlsamen2, c2, "HE") {
            for (i = 1; i <= n; i = i + 1) {
                sumr = 0;
                sumri = 0;
                for (j = 1; j <= n; j = j + 1) {
                    sumr += s[i - 1] * abs1(a[(i - 1) + (j - 1) * lda]) * s[j - 1];
                    sumri += abs1(invhilb[(i - 1) + (j - 1) * ldinvhilb]) / (s[j - 1] * s[i - 1]);
                }
                rnorm = max(rnorm, sumr);
                rinorm = max(rinorm, sumri);
            }
        } else if (Mlsamen2, c2, "GE" || Mlsamen2, c2, "GB") {
            for (i = 1; i <= n; i = i + 1) {
                sumr = 0;
                sumri = 0;
                for (j = 1; j <= n; j = j + 1) {
                    sumr += r[i - 1] * abs1(a[(i - 1) + (j - 1) * lda]) * c[j - 1];
                    sumri += abs1(invhilb[(i - 1) + (j - 1) * ldinvhilb]) / (r[j - 1] * c[i - 1]);
                }
                rnorm = max(rnorm, sumr);
                rinorm = max(rinorm, sumri);
            }
        }
        //
        rnorm = rnorm / abs1(a[(1 - 1)]);
        rcond = 1.0 / (rnorm * rinorm);
        //
        //        Calculating the R for normwise rcond.
        for (i = 1; i <= n; i = i + 1) {
            rinv[i - 1] = 0.0;
        }
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                rinv[i - 1] += abs1(a[(i - 1) + (j - 1) * lda]);
            }
        }
        //
        //        Calculating the Normwise rcond.
        rinorm = 0.0;
        for (i = 1; i <= n; i = i + 1) {
            sumri = 0.0;
            for (j = 1; j <= n; j = j + 1) {
                sumri += abs1(invhilb[(i - 1) + (j - 1) * ldinvhilb] * rinv[j - 1]);
            }
            rinorm = max(rinorm, sumri);
        }
        //
        //        invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
        //        by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
        ncond = abs1(a[(1 - 1)]) / rinorm;
        //
        condthresh = m * eps;
        errthresh = m * eps;
        //
        for (k = 1; k <= nrhs; k = k + 1) {
            normt = 0.0;
            normdif = 0.0;
            cwise_err = 0.0;
            for (i = 1; i <= n; i = i + 1) {
                normt = max(abs1(invhilb[(i - 1) + (k - 1) * ldinvhilb]), normt);
                normdif = max(abs1(x[(i - 1) + (k - 1) * ldx] - invhilb[(i - 1) + (k - 1) * ldinvhilb]), normdif);
                if (invhilb[(i - 1) + (k - 1) * ldinvhilb] != 0.0) {
                    cwise_err = max(abs1(x[(i - 1) + (k - 1) * ldx] - invhilb[(i - 1) + (k - 1) * ldinvhilb]) / abs1(invhilb[(i - 1) + (k - 1) * ldinvhilb]), cwise_err);
                } else if (x[(i - 1) + (k - 1) * ldx] != 0.0) {
                    cwise_err = Rlamch("OVERFLOW");
                }
            }
            if (normt != 0.0) {
                nwise_err = normdif / normt;
            } else if (normdif != 0.0) {
                nwise_err = Rlamch("OVERFLOW");
            } else {
                nwise_err = 0.0;
            }
            //
            for (i = 1; i <= n; i = i + 1) {
                rinv[i - 1] = 0.0;
            }
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    rinv[i - 1] += abs1(a[(i - 1) + (j - 1) * lda] * invhilb[(j - 1) + (k - 1) * ldinvhilb]);
                }
            }
            rinorm = 0.0;
            for (i = 1; i <= n; i = i + 1) {
                sumri = 0.0;
                for (j = 1; j <= n; j = j + 1) {
                    sumri += abs1(invhilb[(i - 1) + (j - 1) * ldinvhilb] * rinv[j - 1] / invhilb[(i - 1) + (k - 1) * ldinvhilb]);
                }
                rinorm = max(rinorm, sumri);
            }
            //        invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
            //        by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
            ccond = abs1(a[(1 - 1)]) / rinorm;
            //
            //        Forward error bound tests
            nwise_bnd = errbnd_n[(k + (bnd_i - 1) * nrhs) - 1];
            cwise_bnd = errbnd_c[(k + (bnd_i - 1) * nrhs) - 1];
            nwise_rcond = errbnd_n[(k + (cond_i - 1) * nrhs) - 1];
            cwise_rcond = errbnd_c[(k + (cond_i - 1) * nrhs) - 1];
            //            write (*,*) 'nwise : ', n, k, ncond, nwise_rcond,
            //     $           condthresh, ncond.ge.condthresh
            //            write (*,*) 'nwise2: ', k, nwise_bnd, nwise_err, errthresh
            if (ncond >= condthresh) {
                nguar = "YES";
                if (nwise_bnd > errthresh) {
                    tstrat[1 - 1] = 1 / (2.0e+0 * eps);
                } else {
                    if (nwise_bnd != 0.0) {
                        tstrat[1 - 1] = nwise_err / nwise_bnd;
                    } else if (nwise_err != 0.0) {
                        tstrat[1 - 1] = 1 / (16.0f * eps);
                    } else {
                        tstrat[1 - 1] = 0.0;
                    }
                    if (tstrat[1 - 1] > 1.0) {
                        tstrat[1 - 1] = 1 / (4.0e+0 * eps);
                    }
                }
            } else {
                nguar = "NO";
                if (nwise_bnd < 1.0) {
                    tstrat[1 - 1] = 1 / (8.0e+0 * eps);
                } else {
                    tstrat[1 - 1] = 1.0;
                }
            }
            //            write (*,*) 'cwise : ', n, k, ccond, cwise_rcond,
            //     $           condthresh, ccond.ge.condthresh
            //            write (*,*) 'cwise2: ', k, cwise_bnd, cwise_err, errthresh
            if (ccond >= condthresh) {
                cguar = "YES";
                if (cwise_bnd > errthresh) {
                    tstrat[2 - 1] = 1 / (2.0e+0 * eps);
                } else {
                    if (cwise_bnd != 0.0) {
                        tstrat[2 - 1] = cwise_err / cwise_bnd;
                    } else if (cwise_err != 0.0) {
                        tstrat[2 - 1] = 1 / (16.0e+0 * eps);
                    } else {
                        tstrat[2 - 1] = 0.0;
                    }
                    if (tstrat[2 - 1] > 1.0) {
                        tstrat[2 - 1] = 1 / (4.0e+0 * eps);
                    }
                }
            } else {
                cguar = "NO";
                if (cwise_bnd < 1.0) {
                    tstrat[2 - 1] = 1 / (8.0e+0 * eps);
                } else {
                    tstrat[2 - 1] = 1.0;
                }
            }
            //
            //     Backwards error test
            tstrat[3 - 1] = berr[k - 1] / eps;
            //
            //     Condition number tests
            tstrat[4 - 1] = rcond / orcond;
            if (rcond >= condthresh && tstrat[4 - 1] < 1.0) {
                tstrat[4 - 1] = 1.0 / tstrat[4 - 1];
            }
            //
            tstrat[5 - 1] = ncond / nwise_rcond;
            if (ncond >= condthresh && tstrat[5 - 1] < 1.0) {
                tstrat[5 - 1] = 1.0 / tstrat[5 - 1];
            }
            //
            tstrat[6 - 1] = ccond / nwise_rcond;
            if (ccond >= condthresh && tstrat[6 - 1] < 1.0) {
                tstrat[6 - 1] = 1.0 / tstrat[6 - 1];
            }
            //
            for (i = 1; i <= ntests; i = i + 1) {
                if (tstrat[i - 1] > thresh) {
                    if (!printed_guide) {
                        write(6, star);
                        write(6, "(3x,i2,': Normwise guaranteed forward error',/,5x,"
                                 "'Guaranteed case: if norm ( abs( Xc - Xt )',"
                                 "' / norm ( Xt ) .LE. ERRBND( *, nwise_i, bnd_i ), then',/,5x,"
                                 "'ERRBND( *, nwise_i, bnd_i ) .LE. MAX(SQRT(N), 10) * EPS')"),
                            1;
                        write(6, "(3x,i2,': Componentwise guaranteed forward error')"), 2;
                        write(6, "(3x,i2,': Backwards error')"), 3;
                        write(6, "(3x,i2,': Reciprocal condition number')"), 4;
                        write(6, "(3x,i2,': Reciprocal normwise condition number')"), 5;
                        write(6, "(3x,i2,': Raw normwise error estimate')"), 6;
                        write(6, "(3x,i2,': Reciprocal componentwise condition number')"), 7;
                        write(6, "(3x,i2,': Raw componentwise error estimate')"), 8;
                        write(6, star);
                        printed_guide = true;
                    }
                    write(6, "(' Z',a2,'SVXX: N =',i2,', RHS = ',i2,', NWISE GUAR. = ',a,"
                             "', CWISE GUAR. = ',a,' test(',i1,') =',g12.5)"),
                        c2, n, k, nguar, cguar, i, tstrat(i);
                    nfail++;
                }
            }
        }
        //
        //$$$         WRITE(*,*)
        //$$$         WRITE(*,*) 'Normwise Error Bounds'
        //$$$         WRITE(*,*) 'Guaranteed error bound: ',ERRBND(NRHS,nwise_i,bnd_i)
        //$$$         WRITE(*,*) 'Reciprocal condition number: ',ERRBND(NRHS,nwise_i,cond_i)
        //$$$         WRITE(*,*) 'Raw error estimate: ',ERRBND(NRHS,nwise_i,rawbnd_i)
        //$$$         WRITE(*,*)
        //$$$         WRITE(*,*) 'Componentwise Error Bounds'
        //$$$         WRITE(*,*) 'Guaranteed error bound: ',ERRBND(NRHS,cwise_i,bnd_i)
        //$$$         WRITE(*,*) 'Reciprocal condition number: ',ERRBND(NRHS,cwise_i,cond_i)
        //$$$         WRITE(*,*) 'Raw error estimate: ',ERRBND(NRHS,cwise_i,rawbnd_i)
        //$$$         prINTEGER *, 'Info: ', info
        //$$$         WRITE(*,*)
        //         WRITE(*,*) 'TSTRAT: ',TSTRAT
        //
    }
    //
    write(6, star);
    if (nfail > 0) {
        write(6, "(' Z',a2,'SVXX: ',i6,' out of ',i6,"
                 "' tests failed to pass the threshold')"),
            c2, nfail, ntests *n + n_aux_tests;
    } else {
        write(6, "(' Z',a2,'SVXX passed the tests of error bounds')"), c2;
    }
    //     Test ratios.
    //
}

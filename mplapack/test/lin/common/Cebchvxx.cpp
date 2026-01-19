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

// Derived from LAPACK routine ZEBCHVXX.
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

void Cebchvxx(REAL const thresh, fem::str_cref path) {
    common cmn;
    common_write write(cmn);
    INTEGER ldabcopy = (nmax - 1) + (nmax - 1) + 1;
    COMPLEX zdum = 0.0;
    fem::str<1> fact = 'E';
    fem::str<1> uplo = 'U';
    fem::str<1> trans = 'N';
    fem::str<1> equed = 'N';
    REAL eps = Rlamch("Epsilon");
    INTEGER nfail = 0;
    INTEGER n_aux_tests = 0;
    const INTEGER nmax = 10;
    INTEGER lda = nmax;
    INTEGER ldab = (nmax - 1) + (nmax - 1) + 1;
    INTEGER ldafb = 2 * (nmax - 1) + (nmax - 1) + 1;
    fem::str<2> c2 = path(2, 3);
    bool printed_guide = false;
    INTEGER n = 0;
    const INTEGER nparams = 2;
    REAL params[nparams];
    INTEGER kl = 0;
    INTEGER ku = 0;
    INTEGER nrhs = 0;
    REAL m = 0.0;
    COMPLEX a[nmax * nmax];
    COMPLEX invhilb[nmax * nmax];
    COMPLEX b[nmax * nmax];
    COMPLEX work[nmax * 3 * 5];
    INTEGER info = 0;
    COMPLEX acopy[nmax * nmax];
    INTEGER j = 0;
    INTEGER i = 0;
    COMPLEX ab[(nmax - 1) + (nmax - 1) + 1 * nmax];
    COMPLEX abcopy[(nmax - 1) + (nmax - 1) + 1 * nmax];
    COMPLEX af[nmax * nmax];
    INTEGER ipiv[nmax];
    REAL s[nmax];
    COMPLEX x[nmax * nmax];
    REAL orcond = 0.0;
    REAL rpvgrw = 0.0;
    REAL berr[nmax];
    const INTEGER nerrbnd = 3;
    REAL errbnd_n[nmax * 3];
    REAL errbnd_c[nmax * 3];
    REAL rwork[3 * nmax];
    COMPLEX afb[2 * (nmax - 1) + (nmax - 1) + 1 * nmax];
    REAL r[nmax];
    REAL c[nmax];
    REAL rcond = 0.0;
    REAL diff[nmax * nmax];
    REAL rnorm = 0.0;
    REAL rinorm = 0.0;
    REAL sumr = 0.0;
    REAL sumri = 0.0;
    REAL rinv[nmax];
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
    fem::str<1> nguar;
    const INTEGER ntests = 6;
    REAL tstrat[ntests];
    fem::str<1> cguar;
    for (n = 1; n <= nmax; n = n + 1) {
        params[1 - 1] = -1.0;
        params[2 - 1] = -1.0;
        kl = n - 1;
        ku = n - 1;
        nrhs = n;
        m = max(sqrt(castREAL(n)), 10.0);
        Clahilb(n, n, a, lda, invhilb, lda, b, lda, work, info, path);
        Clacpy("ALL", n, n, a, nmax, acopy, nmax);
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= kl + ku + 1; i = i + 1) {
                ab[(i - 1) + (j - 1) * ldab] = COMPLEX(0.0, 0.0);
            }
        }
        for (j = 1; j <= n; j = j + 1) {
            for (i = max((INTEGER)1, j - ku); i <= min(n, j + kl); i = i + 1) {
                ab[((ku + 1 + i - j) - 1) + (j - 1) * ldab] = a[(i - 1) + (j - 1) * nmax];
            }
        }
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= kl + ku + 1; i = i + 1) {
                abcopy[(i - 1) + (j - 1) * ldabcopy] = COMPLEX(0.0, 0.0);
            }
        }
        Clacpy("ALL", kl + ku + 1, n, ab, ldab, abcopy, ldab);
        if (Mlsamen(2, c2.elems, "SY")) {
            Csysvxx(fact, uplo, n, nrhs, acopy, lda, af, lda, ipiv, equed, s, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else if (Mlsamen(2, c2.elems, "PO")) {
            Cposvxx(fact, uplo, n, nrhs, acopy, lda, af, lda, equed, s, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else if (Mlsamen(2, c2.elems, "HE")) {
            Chesvxx(fact, uplo, n, nrhs, acopy, lda, af, lda, ipiv, equed, s, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else if (Mlsamen(2, c2.elems, "GB")) {
            Cgbsvxx(fact, trans, n, kl, ku, nrhs, abcopy, ldab, afb, ldafb, ipiv, equed, r, c, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        } else {
            Cgesvxx(fact, trans, n, nrhs, acopy, lda, af, lda, ipiv, equed, r, c, b, lda, x, lda, orcond, rpvgrw, berr, nerrbnd, errbnd_n, errbnd_c, nparams, params, work, rwork, info);
        }
        n_aux_tests++;
        if (orcond < eps) {
        } else {
            if (info > 0 && info <= n + 1) {
                nfail++;
                write(6, "(' Z',a2,'SVXX: N =',i2,', INFO = ',i3,', ORCOND = ',g12.5,"
                         "', real RCOND = ',g12.5)"),
                    c2, n, info, orcond, rcond;
            }
        }
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= nrhs; j = j + 1) {
                diff[(i - 1) + (j - 1) * nmax] = (x[(i - 1) + (j - 1) * nmax] - invhilb[(i - 1) + (j - 1) * nmax]).real();
            }
        }
        rnorm = 0.0;
        rinorm = 0.0;
        if (Mlsamen(2, c2.elems, "PO") || Mlsamen(2, c2.elems, "SY") || Mlsamen(2, c2.elems, "HE")) {
            for (i = 1; i <= n; i = i + 1) {
                sumr = 0.0;
                sumri = 0.0;
                for (j = 1; j <= n; j = j + 1) {
                    sumr += s[i - 1] * cabs1(a[(i - 1) + (j - 1) * nmax]) * s[j - 1];
                    sumri += cabs1(invhilb[(i - 1) + (j - 1) * nmax]) / (s[j - 1] * s[i - 1]);
                }
                rnorm = max(rnorm, sumr);
                rinorm = max(rinorm, sumri);
            }
        } else if (Mlsamen(2, c2.elems, "GE") || Mlsamen(2, c2.elems, "GB")) {
            for (i = 1; i <= n; i = i + 1) {
                sumr = 0.0;
                sumri = 0.0;
                for (j = 1; j <= n; j = j + 1) {
                    sumr += r[i - 1] * cabs1(a[(i - 1) + (j - 1) * nmax]) * c[j - 1];
                    sumri += cabs1(invhilb[(i - 1) + (j - 1) * nmax]) / (r[j - 1] * c[i - 1]);
                }
                rnorm = max(rnorm, sumr);
                rinorm = max(rinorm, sumri);
            }
        }
        rnorm = rnorm / cabs1(a[0]);
        rcond = 1.0 / (rnorm * rinorm);
        for (i = 1; i <= n; i = i + 1) {
            rinv[i - 1] = 0.0;
        }
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                rinv[i - 1] += cabs1(a[(i - 1) + (j - 1) * nmax]);
            }
        }
        rinorm = 0.0;
        for (i = 1; i <= n; i = i + 1) {
            sumri = 0.0;
            for (j = 1; j <= n; j = j + 1) {
                sumri += cabs1(invhilb[(i - 1) + (j - 1) * nmax] * rinv[j - 1]);
            }
            rinorm = max(rinorm, sumri);
        }
        ncond = cabs1(a[0]) / rinorm;
        condthresh = m * eps;
        errthresh = m * eps;
        for (k = 1; k <= nrhs; k = k + 1) {
            normt = 0.0;
            normdif = 0.0;
            cwise_err = 0.0;
            for (i = 1; i <= n; i = i + 1) {
                normt = max(cabs1(invhilb[(i - 1) + (k - 1) * nmax]), normt);
                normdif = max(cabs1(x[(i - 1) + (k - 1) * nmax] - invhilb[(i - 1) + (k - 1) * nmax]), normdif);
                if (invhilb[(i - 1) + (k - 1) * nmax] != 0.0) {
                    cwise_err = max(cabs1(x[(i - 1) + (k - 1) * nmax] - invhilb[(i - 1) + (k - 1) * nmax]) / cabs1(invhilb[(i - 1) + (k - 1) * nmax]), cwise_err);
                } else if (x[(i - 1) + (k - 1) * nmax] != 0.0) {
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
            for (i = 1; i <= n; i = i + 1) {
                rinv[i - 1] = 0.0;
            }
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    rinv[i - 1] += cabs1(a[(i - 1) + (j - 1) * nmax] * invhilb[(j - 1) + (k - 1) * nmax]);
                }
            }
            rinorm = 0.0;
            for (i = 1; i <= n; i = i + 1) {
                sumri = 0.0;
                for (j = 1; j <= n; j = j + 1) {
                    sumri += cabs1(invhilb[(i - 1) + (j - 1) * nmax] * rinv[j - 1] / invhilb[(i - 1) + (k - 1) * nmax]);
                }
                rinorm = max(rinorm, sumri);
            }
            ccond = cabs1(a[0]) / rinorm;
            nwise_bnd = errbnd_n[(k + (bnd_i - 1) * nrhs) - 1];
            cwise_bnd = errbnd_c[(k + (bnd_i - 1) * nrhs) - 1];
            nwise_rcond = errbnd_n[(k + (cond_i - 1) * nrhs) - 1];
            cwise_rcond = errbnd_c[(k + (cond_i - 1) * nrhs) - 1];
            if (ncond >= condthresh) {
                nguar = "YES";
                if (nwise_bnd > errthresh) {
                    tstrat[1 - 1] = 1 / (2.0 * eps);
                } else {
                    if (nwise_bnd != 0.0) {
                        tstrat[1 - 1] = nwise_err / nwise_bnd;
                    } else if (nwise_err != 0.0) {
                        tstrat[1 - 1] = 1 / (16.0 * eps);
                    } else {
                        tstrat[1 - 1] = 0.0;
                    }
                    if (tstrat[1 - 1] > 1.0) {
                        tstrat[1 - 1] = 1 / (4.0 * eps);
                    }
                }
            } else {
                nguar = "NO";
                if (nwise_bnd < 1.0) {
                    tstrat[1 - 1] = 1 / (8.0 * eps);
                } else {
                    tstrat[1 - 1] = 1.0;
                }
            }
            if (ccond >= condthresh) {
                cguar = "YES";
                if (cwise_bnd > errthresh) {
                    tstrat[2 - 1] = 1 / (2.0 * eps);
                } else {
                    if (cwise_bnd != 0.0) {
                        tstrat[2 - 1] = cwise_err / cwise_bnd;
                    } else if (cwise_err != 0.0) {
                        tstrat[2 - 1] = 1 / (16.0 * eps);
                    } else {
                        tstrat[2 - 1] = 0.0;
                    }
                    if (tstrat[2 - 1] > 1.0) {
                        tstrat[2 - 1] = 1 / (4.0 * eps);
                    }
                }
            } else {
                cguar = "NO";
                if (cwise_bnd < 1.0) {
                    tstrat[2 - 1] = 1 / (8.0 * eps);
                } else {
                    tstrat[2 - 1] = 1.0;
                }
            }
            tstrat[3 - 1] = berr[k - 1] / eps;
            tstrat[4 - 1] = rcond / orcond;
            if (rcond >= condthresh && tstrat[4 - 1] < 1.0) {
                tstrat[4 - 1] = 1.0 / tstrat[4 - 1];
            }
            tstrat[5 - 1] = ncond / nwise_rcond;
            if (ncond >= condthresh && tstrat[5 - 1] < 1.0) {
                tstrat[5 - 1] = 1.0 / tstrat[5 - 1];
            }
            tstrat[6 - 1] = ccond / nwise_rcond;
            if (ccond >= condthresh && tstrat[6 - 1] < 1.0) {
                tstrat[6 - 1] = 1.0 / tstrat[6 - 1];
            }
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
                        c2, n, k, nguar, cguar, i, tstrat[i - 1];
                    nfail++;
                }
            }
        }
    }
    write(6, star);
    if (nfail > 0) {
        write(6, "(' Z',a2,'SVXX: ',i6,' out of ',i6,"
                 "' tests failed to pass the threshold')"),
            c2, nfail, ntests *n + n_aux_tests;
    } else {
        write(6, "(' Z',a2,'SVXX passed the tests of error bounds')"), c2;
    }
    // Test ratios.
    //
}

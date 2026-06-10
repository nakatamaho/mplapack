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

// Derived from LAPACK routine ZCKCSD.
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

void Clacsg(INTEGER const m, INTEGER const p, INTEGER const q, REAL *theta, INTEGER (&iseed)[4], COMPLEX *x, INTEGER const ldx, COMPLEX *work) {
    //
    INTEGER r = min(p, m - p, q, m - q);
    //
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    Claset("Full", m, m, zero, zero, x, ldx);
    //
    INTEGER i = 0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    for (i = 1; i <= min(p, q) - r; i = i + 1) {
        x[(i - 1) + (i - 1) * ldx] = one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[((min(p, q) - r + i) - 1) + ((min(p, q) - r + i) - 1) * ldx] = COMPLEX(cos(theta[i - 1]), 0.0);
    }
    for (i = 1; i <= min(p, m - q) - r; i = i + 1) {
        x[((p - i + 1) - 1) + ((m - i + 1) - 1) * ldx] = -one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[((p - (min(p, m - q) - r) + 1 - i) - 1) + ((m - (min(p, m - q) - r) + 1 - i) - 1) * ldx] = COMPLEX(-sin(theta[(r - i + 1) - 1]), 0.0);
    }
    for (i = 1; i <= min(m - p, q) - r; i = i + 1) {
        x[((m - i + 1) - 1) + ((q - i + 1) - 1) * ldx] = one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[((m - (min(m - p, q) - r) + 1 - i) - 1) + ((q - (min(m - p, q) - r) + 1 - i) - 1) * ldx] = COMPLEX(sin(theta[(r - i + 1) - 1]), 0.0);
    }
    for (i = 1; i <= min(m - p, m - q) - r; i = i + 1) {
        x[((p + i) - 1) + ((q + i) - 1) * ldx] = one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x[((p + (min(m - p, m - q) - r) + i) - 1) + ((q + (min(m - p, m - q) - r) + i) - 1) * ldx] = COMPLEX(cos(theta[i - 1]), 0.0);
    }
    INTEGER info = 0;
    Claror("Left", "No init", p, m, x, ldx, iseed, work, info);
    Claror("Left", "No init", m - p, m, &x[((p + 1) - 1)], ldx, iseed, work, info);
    Claror("Right", "No init", m, q, x, ldx, iseed, work, info);
    Claror("Right", "No init", m, m - q, &x[((q + 1) - 1) * ldx], ldx, iseed, work, info);
    //
}

void Cckcsd(INTEGER const nm, INTEGER *mval, INTEGER *pval, INTEGER *qval, INTEGER const nmats, INTEGER (&iseed)[4], REAL const thresh, INTEGER const mmax, COMPLEX *x, COMPLEX *xf, COMPLEX *u1, COMPLEX *u2, COMPLEX *v1t, COMPLEX *v2t, REAL *theta, INTEGER *iwork, COMPLEX *work, REAL *rwork, INTEGER const nin, INTEGER const nout, INTEGER &info) {
    common cmn;
    common_write write(cmn);
    fem::str<3> path;
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    bool firstt = false;
    const INTEGER ntypes = 4;
    bool dotype[ntypes];
    INTEGER ldx = 0;
    INTEGER ldu1 = 0;
    INTEGER ldu2 = 0;
    INTEGER ldv1t = 0;
    INTEGER ldv2t = 0;
    INTEGER lwork = 0;
    INTEGER im = 0;
    INTEGER m = 0;
    INTEGER p = 0;
    INTEGER q = 0;
    INTEGER imat = 0;
    INTEGER iinfo = 0;
    INTEGER r = 0;
    INTEGER i = 0;
    const REAL piover2 = 1.5707963267948966192313216916397514421;
    INTEGER j = 0;
    const REAL orth = 0.000000000001;
    const REAL ten = 10.0;
    const REAL gapdigit = 18.0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    const REAL realzero = 0.0;
    const REAL realone = 1.0;
    INTEGER nt = 0;
    const INTEGER ntests = 15;
    REAL result[ntests];
    //
    static const char *format_9999 = "(' Claror in Cckcsd: M = ',i5,', INFO = ',i15)";
    static const char *format_9998 = "(' M=',i4,' P=',i4,', Q=',i4,', type ',i2,', test ',i2,', ratio=',g13.6)";
    //
    // Initialize constants and the random number seed.
    //
    path(1, 3) = "CSD";
    info = 0;
    nrun = 0;
    nfail = 0;
    firstt = true;
    Alareq(path, nmats, dotype, ntypes, nin, nout);
    ldx = mmax;
    ldu1 = mmax;
    ldu2 = mmax;
    ldv1t = mmax;
    ldv2t = mmax;
    lwork = mmax * mmax;
    //
    // Do for each value of M in MVAL.
    //
    for (im = 1; im <= nm; im = im + 1) {
        m = mval[im - 1];
        p = pval[im - 1];
        q = qval[im - 1];
        //
        for (imat = 1; imat <= ntypes; imat = imat + 1) {
            //
            // Do the tests only if DOTYPE( IMAT ) is true.
            //
            if (!dotype[imat - 1]) {
                goto statement_20;
            }
            //
            // Generate X
            //
            if (imat == 1) {
                Claror("L", "I", m, m, x, ldx, iseed, work, iinfo);
                if (m != 0 && iinfo != 0) {
                    write(nout, format_9999), m, iinfo;
                    info = abs(iinfo);
                    goto statement_20;
                }
            } else if (imat == 2) {
                r = min(p, m - p, q, m - q);
                for (i = 1; i <= r; i = i + 1) {
                    theta[i - 1] = piover2 * Rlarnd(1, iseed);
                }
                Clacsg(m, p, q, theta, iseed, x, ldx, work);
                for (i = 1; i <= m; i = i + 1) {
                    for (j = 1; j <= m; j = j + 1) {
                        x[(i + (j - 1) * ldx) - 1] += orth * Rlarnd(2, iseed);
                    }
                }
            } else if (imat == 3) {
                r = min(p, m - p, q, m - q);
                for (i = 1; i <= r + 1; i = i + 1) {
                    theta[i - 1] = pow(ten, (-Rlarnd(1, iseed) * gapdigit));
                }
                for (i = 2; i <= r + 1; i = i + 1) {
                    theta[i - 1] += theta[(i - 1) - 1];
                }
                for (i = 1; i <= r; i = i + 1) {
                    theta[i - 1] = piover2 * theta[i - 1] / theta[(r + 1) - 1];
                }
                Clacsg(m, p, q, theta, iseed, x, ldx, work);
            } else {
                Claset("F", m, m, zero, one, x, ldx);
                for (i = 1; i <= m; i = i + 1) {
                    j = castINTEGER(Rlaran(iseed) * m) + 1;
                    if (j != i) {
                        CRrot(m, &x[(1 + (i - 1) * ldx) - 1], 1, &x[(1 + (j - 1) * ldx) - 1], 1, realzero, realone);
                    }
                }
            }
            //
            nt = 15;
            //
            Ccsdts(m, p, q, x, xf, ldx, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, theta, iwork, work, lwork, rwork, result);
            //
            // Print information about the tests that did not
            // pass the threshold.
            //
            for (i = 1; i <= nt; i = i + 1) {
                if (result[i - 1] >= thresh) {
                    if (nfail == 0 && firstt) {
                        firstt = false;
                        Alahdg(nout, path);
                    }
                    write(nout, format_9998), m, p, q, imat, i, result[i - 1];
                    nfail++;
                }
            }
            nrun += nt;
        statement_20:;
        }
    }
    //
    // Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, 0);
    //
    // End of Cckcsd
    //
}

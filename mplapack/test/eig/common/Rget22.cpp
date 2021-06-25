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

void Rget22(const char *transa, const char *transe, const char *transw, INTEGER const n, REAL *a, INTEGER const lda, REAL *e, INTEGER const lde, REAL *wr, REAL *wi, REAL *work, REAL *result) {
    //
    //     Initialize RESULT (in case N=0)
    //
    const REAL zero = 0.0;
    result[1 - 1] = zero;
    result[2 - 1] = zero;
    if (n <= 0) {
        return;
    }
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Precision");
    //
    INTEGER itrnse = 0;
    INTEGER ince = 1;
    char norma;
    char norme;
    norma = 'O';
    norme = 'O';
    //
    if (Mlsame(transa, "T") || Mlsame(transa, "C")) {
        norma = 'I';
    }
    if (Mlsame(transe, "T") || Mlsame(transe, "C")) {
        norme = 'I';
        itrnse = 1;
        ince = lde;
    }
    //
    //     Check normalization of E
    //
    const REAL one = 1.0;
    REAL enrmin = one / ulp;
    REAL enrmax = zero;
    INTEGER ipair = 0;
    INTEGER jvec = 0;
    REAL temp1 = 0.0;
    INTEGER j = 0;
    if (itrnse == 0) {
        //
        //        Eigenvectors are column vectors.
        //
        ipair = 0;
        for (jvec = 1; jvec <= n; jvec = jvec + 1) {
            temp1 = zero;
            if (ipair == 0 && jvec < n && wi[jvec - 1] != zero) {
                ipair = 1;
            }
            if (ipair == 1) {
                //
                //              Complex eigenvector
                //
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max(temp1, abs(e[(j - 1) + (jvec - 1) * lde]) + abs(e[(j - 1) + ((jvec + 1) - 1) * lde]));
                }
                enrmin = min(enrmin, temp1);
                enrmax = max(enrmax, temp1);
                ipair = 2;
            } else if (ipair == 2) {
                ipair = 0;
            } else {
                //
                //              Real eigenvector
                //
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max(temp1, abs(e[(j - 1) + (jvec - 1) * lde]));
                }
                enrmin = min(enrmin, temp1);
                enrmax = max(enrmax, temp1);
                ipair = 0;
            }
        }
        //
    } else {
        //
        //        Eigenvectors are row vectors.
        //
        for (jvec = 1; jvec <= n; jvec = jvec + 1) {
            work[jvec - 1] = zero;
        }
        //
        for (j = 1; j <= n; j = j + 1) {
            ipair = 0;
            for (jvec = 1; jvec <= n; jvec = jvec + 1) {
                if (ipair == 0 && jvec < n && wi[jvec - 1] != zero) {
                    ipair = 1;
                }
                if (ipair == 1) {
                    work[jvec - 1] = max(work[jvec - 1], abs(e[(j - 1) + (jvec - 1) * lde]) + abs(e[(j - 1) + ((jvec + 1) - 1) * lde]));
                    work[(jvec + 1) - 1] = work[jvec - 1];
                } else if (ipair == 2) {
                    ipair = 0;
                } else {
                    work[jvec - 1] = max(work[jvec - 1], abs(e[(j - 1) + (jvec - 1) * lde]));
                    ipair = 0;
                }
            }
        }
        //
        for (jvec = 1; jvec <= n; jvec = jvec + 1) {
            enrmin = min(enrmin, work[jvec - 1]);
            enrmax = max(enrmax, work[jvec - 1]);
        }
    }
    //
    //     Norm of A:
    //
    REAL anorm = max(Rlange(&norma, n, n, a, lda, work), unfl);
    //
    //     Norm of E:
    //
    REAL enorm = max(Rlange(&norme, n, n, e, lde, work), ulp);
    //
    //     Norm of error:
    //
    //     Error =  AE - EW
    //
    Rlaset("Full", n, n, zero, zero, work, n);
    //
    ipair = 0;
    INTEGER ierow = 1;
    INTEGER iecol = 1;
    //
    INTEGER jcol = 0;
    REAL wmat[2 * 2];
    INTEGER ldwmat = 2;
    for (jcol = 1; jcol <= n; jcol = jcol + 1) {
        if (itrnse == 1) {
            ierow = jcol;
        } else {
            iecol = jcol;
        }
        //
        if (ipair == 0 && wi[jcol - 1] != zero) {
            ipair = 1;
        }
        //
        if (ipair == 1) {
            wmat[(1 - 1)] = wr[jcol - 1];
            wmat[(2 - 1)] = -wi[jcol - 1];
            wmat[(2 - 1) * ldwmat] = wi[jcol - 1];
            wmat[(2 - 1) + (2 - 1) * ldwmat] = wr[jcol - 1];
            Rgemm(transe, transw, n, 2, 2, one, &e[(ierow - 1) + (iecol - 1) * lde], lde, wmat, 2, zero, &work[(n * (jcol - 1) + 1) - 1], n);
            ipair = 2;
        } else if (ipair == 2) {
            ipair = 0;
            //
        } else {
            //
            Raxpy(n, wr[jcol - 1], &e[(ierow - 1) + (iecol - 1) * lde], ince, &work[(n * (jcol - 1) + 1) - 1], 1);
            ipair = 0;
        }
        //
    }
    //
    Rgemm(transa, transe, n, n, n, one, a, lda, e, lde, -one, work, n);
    //
    REAL errnrm = Rlange("One", n, n, work, n, &work[(n * n + 1) - 1]) / enorm;
    //
    //     Compute RESULT(1) (avoiding under/overflow)
    //
    if (anorm > errnrm) {
        result[1 - 1] = (errnrm / anorm) / ulp;
    } else {
        if (anorm < one) {
            result[1 - 1] = (min(errnrm, anorm) / anorm) / ulp;
        } else {
            result[1 - 1] = min(errnrm / anorm, one) / ulp;
        }
    }
    //
    //     Compute RESULT(2) : the normalization error in E.
    //
    result[2 - 1] = max(abs(enrmax - one), abs(enrmin - one)) / (castREAL(n) * ulp);
    //
    //     End of Rget22
    //
}

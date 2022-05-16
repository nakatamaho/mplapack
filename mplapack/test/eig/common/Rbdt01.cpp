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

void Rbdt01(INTEGER const m, INTEGER const n, INTEGER const kd, REAL *a, INTEGER const lda, REAL *q, INTEGER const ldq, REAL *d, REAL *e, REAL *pt, INTEGER const ldpt, REAL *work, REAL &resid) {
    //
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    if (m <= 0 || n <= 0) {
        resid = zero;
        return;
    }
    //
    //     Compute A - Q * B * P' one column at a time.
    //
    resid = zero;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL one = 1.0;
    if (kd != 0) {
        //
        //        B is bidiagonal.
        //
        if (kd != 0 && m >= n) {
            //
            //           B is upper bidiagonal and M >= N.
            //
            for (j = 1; j <= n; j = j + 1) {
                Rcopy(m, &a[(j - 1) * lda], 1, work, 1);
                for (i = 1; i <= n - 1; i = i + 1) {
                    work[(m + i) - 1] = d[i - 1] * pt[(i - 1) + (j - 1) * ldpt] + e[i - 1] * pt[((i + 1) - 1) + (j - 1) * ldpt];
                }
                work[(m + n) - 1] = d[n - 1] * pt[(n - 1) + (j - 1) * ldpt];
                Rgemv("No transpose", m, n, -one, q, ldq, &work[(m + 1) - 1], 1, one, work, 1);
                resid = max({resid, Rasum(m, work, 1)});
            }
        } else if (kd < 0) {
            //
            //           B is upper bidiagonal and M < N.
            //
            for (j = 1; j <= n; j = j + 1) {
                Rcopy(m, &a[(j - 1) * lda], 1, work, 1);
                for (i = 1; i <= m - 1; i = i + 1) {
                    work[(m + i) - 1] = d[i - 1] * pt[(i - 1) + (j - 1) * ldpt] + e[i - 1] * pt[((i + 1) - 1) + (j - 1) * ldpt];
                }
                work[(m + m) - 1] = d[m - 1] * pt[(m - 1) + (j - 1) * ldpt];
                Rgemv("No transpose", m, m, -one, q, ldq, &work[(m + 1) - 1], 1, one, work, 1);
                resid = max({resid, Rasum(m, work, 1)});
            }
        } else {
            //
            //           B is lower bidiagonal.
            //
            for (j = 1; j <= n; j = j + 1) {
                Rcopy(m, &a[(j - 1) * lda], 1, work, 1);
                work[(m + 1) - 1] = d[1 - 1] * pt[(j - 1) * ldpt];
                for (i = 2; i <= m; i = i + 1) {
                    work[(m + i) - 1] = e[(i - 1) - 1] * pt[((i - 1) - 1) + (j - 1) * ldpt] + d[i - 1] * pt[(i - 1) + (j - 1) * ldpt];
                }
                Rgemv("No transpose", m, m, -one, q, ldq, &work[(m + 1) - 1], 1, one, work, 1);
                resid = max({resid, Rasum(m, work, 1)});
            }
        }
    } else {
        //
        //        B is diagonal.
        //
        if (m >= n) {
            for (j = 1; j <= n; j = j + 1) {
                Rcopy(m, &a[(j - 1) * lda], 1, work, 1);
                for (i = 1; i <= n; i = i + 1) {
                    work[(m + i) - 1] = d[i - 1] * pt[(i - 1) + (j - 1) * ldpt];
                }
                Rgemv("No transpose", m, n, -one, q, ldq, &work[(m + 1) - 1], 1, one, work, 1);
                resid = max({resid, Rasum(m, work, 1)});
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                Rcopy(m, &a[(j - 1) * lda], 1, work, 1);
                for (i = 1; i <= m; i = i + 1) {
                    work[(m + i) - 1] = d[i - 1] * pt[(i - 1) + (j - 1) * ldpt];
                }
                Rgemv("No transpose", m, m, -one, q, ldq, &work[(m + 1) - 1], 1, one, work, 1);
                resid = max({resid, Rasum(m, work, 1)});
            }
        }
    }
    //
    //     Compute norm(A - Q * B * P') / ( n * norm(A) * EPS )
    //
    REAL anorm = Rlange("1", m, n, a, lda, work);
    REAL eps = Rlamch("Precision");
    //
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        if (anorm >= resid) {
            resid = (resid / anorm) / (castREAL(n) * eps);
        } else {
            if (anorm < one) {
                resid = (min(resid, REAL(castREAL(n) * anorm)) / anorm) / (castREAL(n) * eps);
            } else {
                resid = min(REAL(resid / anorm), castREAL(n)) / (castREAL(n) * eps);
            }
        }
    }
    //
    //     End of Rbdt01
    //
}

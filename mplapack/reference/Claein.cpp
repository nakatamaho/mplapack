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

inline REAL abs1(COMPLEX cdum) { return abs(cdum.real()) + abs(cdum.imag()); }

void Claein(bool const rightv, bool const noinit, INTEGER const n, COMPLEX *h, INTEGER const ldh, COMPLEX const w, COMPLEX *v, COMPLEX *b, INTEGER const ldb, REAL *rwork, REAL const eps3, REAL const smlnum, INTEGER &info) {
    COMPLEX cdum = 0.0;
    REAL rootn = 0.0;
    const REAL tenth = 1.0e-1;
    REAL growto = 0.0;
    const REAL one = 1.0;
    REAL nrmsml = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    REAL vnorm = 0.0;
    COMPLEX ei = 0.0;
    COMPLEX x = 0.0;
    COMPLEX temp = 0.0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    char trans;
    COMPLEX ej = 0.0;
    char normin;
    INTEGER its = 0;
    REAL scale = 0.0;
    INTEGER ierr = 0;
    REAL rtemp = 0.0;
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    //
    //     GROWTO is the threshold used in the acceptance test for an
    //     eigenvector.
    //
    rootn = sqrt(castREAL(n));
    growto = tenth / rootn;
    nrmsml = max(one, eps3 * rootn) * smlnum;
    //
    //     Form B = H - W*I (except that the subdiagonal elements are not
    //     stored).
    //
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= j - 1; i = i + 1) {
            b[(i - 1) + (j - 1) * ldb] = h[(i - 1) + (j - 1) * ldh];
        }
        b[(j - 1) + (j - 1) * ldb] = h[(j - 1) + (j - 1) * ldh] - w;
    }
    //
    if (noinit) {
        //
        //        Initialize V.
        //
        for (i = 1; i <= n; i = i + 1) {
            v[i - 1] = eps3;
        }
    } else {
        //
        //        Scale supplied initial vector.
        //
        vnorm = RCnrm2(n, v, 1);
        CRscal(n, (eps3 * rootn) / max(vnorm, nrmsml), v, 1);
    }
    //
    if (rightv) {
        //
        //        LU decomposition with partial pivoting of B, replacing zero
        //        pivots by EPS3.
        //
        for (i = 1; i <= n - 1; i = i + 1) {
            ei = h[((i + 1) - 1) + (i - 1) * ldh];
            if (abs1(b[(i - 1) + (i - 1) * ldb]) < abs1(ei)) {
                //
                //              Interchange rows and eliminate.
                //
                x = Cladiv(b[(i - 1) + (i - 1) * ldb], ei);
                b[(i - 1) + (i - 1) * ldb] = ei;
                for (j = i + 1; j <= n; j = j + 1) {
                    temp = b[((i + 1) - 1) + (j - 1) * ldb];
                    b[((i + 1) - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - x * temp;
                    b[(i - 1) + (j - 1) * ldb] = temp;
                }
            } else {
                //
                //              Eliminate without interchange.
                //
                if (b[(i - 1) + (i - 1) * ldb] == zero) {
                    b[(i - 1) + (i - 1) * ldb] = eps3;
                }
                x = Cladiv(ei, b[(i - 1) + (i - 1) * ldb]);
                if (x != zero) {
                    for (j = i + 1; j <= n; j = j + 1) {
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb] - x * b[(i - 1) + (j - 1) * ldb];
                    }
                }
            }
        }
        if (b[(n - 1) + (n - 1) * ldb] == zero) {
            b[(n - 1) + (n - 1) * ldb] = eps3;
        }
        //
        trans = 'N';
        //
    } else {
        //
        //        UL decomposition with partial pivoting of B, replacing zero
        //        pivots by EPS3.
        //
        for (j = n; j >= 2; j = j - 1) {
            ej = h[(j - 1) + ((j - 1) - 1) * ldh];
            if (abs1(b[(j - 1) + (j - 1) * ldb]) < abs1(ej)) {
                //
                //              Interchange columns and eliminate.
                //
                x = Cladiv(b[(j - 1) + (j - 1) * ldb], ej);
                b[(j - 1) + (j - 1) * ldb] = ej;
                for (i = 1; i <= j - 1; i = i + 1) {
                    temp = b[(i - 1) + ((j - 1) - 1) * ldb];
                    b[(i - 1) + ((j - 1) - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - x * temp;
                    b[(i - 1) + (j - 1) * ldb] = temp;
                }
            } else {
                //
                //              Eliminate without interchange.
                //
                if (b[(j - 1) + (j - 1) * ldb] == zero) {
                    b[(j - 1) + (j - 1) * ldb] = eps3;
                }
                x = Cladiv(ej, b[(j - 1) + (j - 1) * ldb]);
                if (x != zero) {
                    for (i = 1; i <= j - 1; i = i + 1) {
                        b[(i - 1) + ((j - 1) - 1) * ldb] = b[(i - 1) + ((j - 1) - 1) * ldb] - x * b[(i - 1) + (j - 1) * ldb];
                    }
                }
            }
        }
        if (b[(1 - 1)] == zero) {
            b[(1 - 1)] = eps3;
        }
        //
        trans = 'C';
        //
    }
    //
    normin = 'N';
    for (its = 1; its <= n; its = its + 1) {
        //
        //        Solve U*x = scale*v for a right eigenvector
        //          or U**H *x = scale*v for a left eigenvector,
        //        overwriting x on v.
        //
        Clatrs("Upper", &trans, "Nonunit", &normin, n, b, ldb, v, scale, rwork, ierr);
        normin = 'Y';
        //
        //        Test for sufficient growth in the norm of v.
        //
        vnorm = RCasum(n, v, 1);
        if (vnorm >= growto * scale) {
            goto statement_120;
        }
        //
        //        Choose new orthogonal starting vector and try again.
        //
        rtemp = eps3 / (rootn + one);
        v[1 - 1] = eps3;
        for (i = 2; i <= n; i = i + 1) {
            v[i - 1] = rtemp;
        }
        v[(n - its + 1) - 1] = v[(n - its + 1) - 1] - eps3 * rootn;
    }
    //
    //     Failure to find eigenvector in N iterations.
    //
    info = 1;
//
statement_120:
    //
    //     Normalize eigenvector.
    //
    i = iCamax(n, v, 1);
    CRscal(n, one / abs1(v[i - 1]), v, 1);
    //
    //     End of Claein
    //
}

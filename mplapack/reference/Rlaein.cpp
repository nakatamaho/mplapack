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

void Rlaein(bool const rightv, bool const noinit, INTEGER const n, REAL *h, INTEGER const ldh, REAL const wr, REAL const wi, REAL *vr, REAL *vi, REAL *b, INTEGER const ldb, REAL *work, REAL const eps3, REAL const smlnum, REAL const bignum, INTEGER &info) {
    REAL rootn = 0.0;
    const REAL tenth = 1.0e-1;
    REAL growto = 0.0;
    const REAL one = 1.0;
    REAL nrmsml = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL vnorm = 0.0;
    REAL ei = 0.0;
    REAL x = 0.0;
    REAL temp = 0.0;
    char trans;
    REAL ej = 0.0;
    char normin;
    INTEGER its = 0;
    REAL scale = 0.0;
    INTEGER ierr = 0;
    REAL norm = 0.0;
    REAL rec = 0.0;
    REAL absbii = 0.0;
    REAL xr = 0.0;
    REAL xi = 0.0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    REAL absbjj = 0.0;
    REAL vmax = 0.0;
    REAL vcrit = 0.0;
    REAL w = 0.0;
    REAL w1 = 0.0;
    REAL y = 0.0;
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
    //     .. Executable Statements ..
    //
    info = 0;
    //
    //     GROWTO is the threshold used in the acceptance test for an
    //     eigenvector.
    //
    rootn = sqrt(n.real());
    growto = tenth / rootn;
    nrmsml = max(one, eps3 * rootn) * smlnum;
    //
    //     Form B = H - (WR,WI)*I (except that the subdiagonal elements and
    //     the imaginary parts of the diagonal elements are not stored).
    //
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= j - 1; i = i + 1) {
            b[(i - 1) + (j - 1) * ldb] = h[(i - 1) + (j - 1) * ldh];
        }
        b[(j - 1) + (j - 1) * ldb] = h[(j - 1) + (j - 1) * ldh] - wr;
    }
    //
    if (wi == zero) {
        //
        //        Real eigenvalue.
        //
        if (noinit) {
            //
            //           Set initial vector.
            //
            for (i = 1; i <= n; i = i + 1) {
                vr[i - 1] = eps3;
            }
        } else {
            //
            //           Scale supplied initial vector.
            //
            vnorm = Rnrm2(n, vr, 1);
            Rscal(n, (eps3 * rootn) / max(vnorm, nrmsml), vr, 1);
        }
        //
        if (rightv) {
            //
            //           LU decomposition with partial pivoting of B, replacing zero
            //           pivots by EPS3.
            //
            for (i = 1; i <= n - 1; i = i + 1) {
                ei = h[((i + 1) - 1) + (i - 1) * ldh];
                if (abs(b[(i - 1) + (i - 1) * ldb]) < abs(ei)) {
                    //
                    //                 Interchange rows and eliminate.
                    //
                    x = b[(i - 1) + (i - 1) * ldb] / ei;
                    b[(i - 1) + (i - 1) * ldb] = ei;
                    for (j = i + 1; j <= n; j = j + 1) {
                        temp = b[((i + 1) - 1) + (j - 1) * ldb];
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - x * temp;
                        b[(i - 1) + (j - 1) * ldb] = temp;
                    }
                } else {
                    //
                    //                 Eliminate without interchange.
                    //
                    if (b[(i - 1) + (i - 1) * ldb] == zero) {
                        b[(i - 1) + (i - 1) * ldb] = eps3;
                    }
                    x = ei / b[(i - 1) + (i - 1) * ldb];
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
            //           UL decomposition with partial pivoting of B, replacing zero
            //           pivots by EPS3.
            //
            for (j = n; j >= 2; j = j - 1) {
                ej = h[(j - 1) + ((j - 1) - 1) * ldh];
                if (abs(b[(j - 1) + (j - 1) * ldb]) < abs(ej)) {
                    //
                    //                 Interchange columns and eliminate.
                    //
                    x = b[(j - 1) + (j - 1) * ldb] / ej;
                    b[(j - 1) + (j - 1) * ldb] = ej;
                    for (i = 1; i <= j - 1; i = i + 1) {
                        temp = b[(i - 1) + ((j - 1) - 1) * ldb];
                        b[(i - 1) + ((j - 1) - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - x * temp;
                        b[(i - 1) + (j - 1) * ldb] = temp;
                    }
                } else {
                    //
                    //                 Eliminate without interchange.
                    //
                    if (b[(j - 1) + (j - 1) * ldb] == zero) {
                        b[(j - 1) + (j - 1) * ldb] = eps3;
                    }
                    x = ej / b[(j - 1) + (j - 1) * ldb];
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
            trans = 'T';
            //
        }
        //
        normin = "N";
        for (its = 1; its <= n; its = its + 1) {
            //
            //           Solve U*x = scale*v for a right eigenvector
            //             or U**T*x = scale*v for a left eigenvector,
            //           overwriting x on v.
            //
            Rlatrs("Upper", trans, "Nonunit", normin, n, b, ldb, vr, scale, work, ierr);
            normin = "Y";
            //
            //           Test for sufficient growth in the norm of v.
            //
            vnorm = Rasum(n, vr, 1);
            if (vnorm >= growto * scale) {
                goto statement_120;
            }
            //
            //           Choose new orthogonal starting vector and try again.
            //
            temp = eps3 / (rootn + one);
            vr[1 - 1] = eps3;
            for (i = 2; i <= n; i = i + 1) {
                vr[i - 1] = temp;
            }
            vr[(n - its + 1) - 1] = vr[(n - its + 1) - 1] - eps3 * rootn;
        }
        //
        //        Failure to find eigenvector in N iterations.
        //
        info = 1;
    //
    statement_120:
        //
        //        Normalize eigenvector.
        //
        i = iRamax(n, vr, 1);
        Rscal(n, one / abs(vr[i - 1]), vr, 1);
    } else {
        //
        //        Complex eigenvalue.
        //
        if (noinit) {
            //
            //           Set initial vector.
            //
            for (i = 1; i <= n; i = i + 1) {
                vr[i - 1] = eps3;
                vi[i - 1] = zero;
            }
        } else {
            //
            //           Scale supplied initial vector.
            //
            norm = Rlapy2(Rnrm2(n, vr, 1), Rnrm2(n, vi, 1));
            rec = (eps3 * rootn) / max(norm, nrmsml);
            Rscal(n, rec, vr, 1);
            Rscal(n, rec, vi, 1);
        }
        //
        if (rightv) {
            //
            //           LU decomposition with partial pivoting of B, replacing zero
            //           pivots by EPS3.
            //
            //           The imaginary part of the (i,j)-th element of U is stored in
            //           B(j+1,i).
            //
            b[(2 - 1)] = -wi;
            for (i = 2; i <= n; i = i + 1) {
                b[((i + 1) - 1)] = zero;
            }
            //
            for (i = 1; i <= n - 1; i = i + 1) {
                absbii = Rlapy2(b[(i - 1) + (i - 1) * ldb], &b[((i + 1) - 1) + (i - 1) * ldb]);
                ei = h[((i + 1) - 1) + (i - 1) * ldh];
                if (absbii < abs(ei)) {
                    //
                    //                 Interchange rows and eliminate.
                    //
                    xr = b[(i - 1) + (i - 1) * ldb] / ei;
                    xi = b[((i + 1) - 1) + (i - 1) * ldb] / ei;
                    b[(i - 1) + (i - 1) * ldb] = ei;
                    b[((i + 1) - 1) + (i - 1) * ldb] = zero;
                    for (j = i + 1; j <= n; j = j + 1) {
                        temp = b[((i + 1) - 1) + (j - 1) * ldb];
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - xr * temp;
                        b[((j + 1) - 1) + ((i + 1) - 1) * ldb] = b[((j + 1) - 1) + (i - 1) * ldb] - xi * temp;
                        b[(i - 1) + (j - 1) * ldb] = temp;
                        b[((j + 1) - 1) + (i - 1) * ldb] = zero;
                    }
                    b[((i + 2) - 1) + (i - 1) * ldb] = -wi;
                    b[((i + 1) - 1) + ((i + 1) - 1) * ldb] = b[((i + 1) - 1) + ((i + 1) - 1) * ldb] - xi * wi;
                    b[((i + 2) - 1) + ((i + 1) - 1) * ldb] += xr * wi;
                } else {
                    //
                    //                 Eliminate without interchanging rows.
                    //
                    if (absbii == zero) {
                        b[(i - 1) + (i - 1) * ldb] = eps3;
                        b[((i + 1) - 1) + (i - 1) * ldb] = zero;
                        absbii = eps3;
                    }
                    ei = (ei / absbii) / absbii;
                    xr = b[(i - 1) + (i - 1) * ldb] * ei;
                    xi = -b[((i + 1) - 1) + (i - 1) * ldb] * ei;
                    for (j = i + 1; j <= n; j = j + 1) {
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb] - xr * b[(i - 1) + (j - 1) * ldb] + xi * b[((j + 1) - 1) + (i - 1) * ldb];
                        b[((j + 1) - 1) + ((i + 1) - 1) * ldb] = -xr * b[((j + 1) - 1) + (i - 1) * ldb] - xi * b[(i - 1) + (j - 1) * ldb];
                    }
                    b[((i + 2) - 1) + ((i + 1) - 1) * ldb] = b[((i + 2) - 1) + ((i + 1) - 1) * ldb] - wi;
                }
                //
                //              Compute 1-norm of offdiagonal elements of i-th row.
                //
                work[i - 1] = Rasum(n - i, &b[(i - 1) + ((i + 1) - 1) * ldb], ldb) + Rasum(n - i, &b[((i + 2) - 1) + (i - 1) * ldb], 1);
            }
            if (b[(n - 1) + (n - 1) * ldb] == zero && b[((n + 1) - 1) + (n - 1) * ldb] == zero) {
                b[(n - 1) + (n - 1) * ldb] = eps3;
            }
            work[n - 1] = zero;
            //
            i1 = n;
            i2 = 1;
            i3 = -1;
        } else {
            //
            //           UL decomposition with partial pivoting of conj(B),
            //           replacing zero pivots by EPS3.
            //
            //           The imaginary part of the (i,j)-th element of U is stored in
            //           B(j+1,i).
            //
            b[((n + 1) - 1) + (n - 1) * ldb] = wi;
            for (j = 1; j <= n - 1; j = j + 1) {
                b[((n + 1) - 1) + (j - 1) * ldb] = zero;
            }
            //
            for (j = n; j >= 2; j = j - 1) {
                ej = h[(j - 1) + ((j - 1) - 1) * ldh];
                absbjj = Rlapy2(b[(j - 1) + (j - 1) * ldb], &b[((j + 1) - 1) + (j - 1) * ldb]);
                if (absbjj < abs(ej)) {
                    //
                    //                 Interchange columns and eliminate
                    //
                    xr = b[(j - 1) + (j - 1) * ldb] / ej;
                    xi = b[((j + 1) - 1) + (j - 1) * ldb] / ej;
                    b[(j - 1) + (j - 1) * ldb] = ej;
                    b[((j + 1) - 1) + (j - 1) * ldb] = zero;
                    for (i = 1; i <= j - 1; i = i + 1) {
                        temp = b[(i - 1) + ((j - 1) - 1) * ldb];
                        b[(i - 1) + ((j - 1) - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - xr * temp;
                        b[(j - 1) + (i - 1) * ldb] = b[((j + 1) - 1) + (i - 1) * ldb] - xi * temp;
                        b[(i - 1) + (j - 1) * ldb] = temp;
                        b[((j + 1) - 1) + (i - 1) * ldb] = zero;
                    }
                    b[((j + 1) - 1) + ((j - 1) - 1) * ldb] = wi;
                    b[((j - 1) - 1) + ((j - 1) - 1) * ldb] += xi * wi;
                    b[(j - 1) + ((j - 1) - 1) * ldb] = b[(j - 1) + ((j - 1) - 1) * ldb] - xr * wi;
                } else {
                    //
                    //                 Eliminate without interchange.
                    //
                    if (absbjj == zero) {
                        b[(j - 1) + (j - 1) * ldb] = eps3;
                        b[((j + 1) - 1) + (j - 1) * ldb] = zero;
                        absbjj = eps3;
                    }
                    ej = (ej / absbjj) / absbjj;
                    xr = b[(j - 1) + (j - 1) * ldb] * ej;
                    xi = -b[((j + 1) - 1) + (j - 1) * ldb] * ej;
                    for (i = 1; i <= j - 1; i = i + 1) {
                        b[(i - 1) + ((j - 1) - 1) * ldb] = b[(i - 1) + ((j - 1) - 1) * ldb] - xr * b[(i - 1) + (j - 1) * ldb] + xi * b[((j + 1) - 1) + (i - 1) * ldb];
                        b[(j - 1) + (i - 1) * ldb] = -xr * b[((j + 1) - 1) + (i - 1) * ldb] - xi * b[(i - 1) + (j - 1) * ldb];
                    }
                    b[(j - 1) + ((j - 1) - 1) * ldb] += wi;
                }
                //
                //              Compute 1-norm of offdiagonal elements of j-th column.
                //
                work[j - 1] = Rasum(j - 1, &b[(j - 1) * ldb], 1) + Rasum(j - 1, &b[((j + 1) - 1)], ldb);
            }
            if (b[(1 - 1)] == zero && b[(2 - 1)] == zero) {
                b[(1 - 1)] = eps3;
            }
            work[1 - 1] = zero;
            //
            i1 = 1;
            i2 = n;
            i3 = 1;
        }
        //
        for (its = 1; its <= n; its = its + 1) {
            scale = one;
            vmax = one;
            vcrit = bignum;
            //
            //           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
            //             or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector,
            //           overwriting (xr,xi) on (vr,vi).
            //
            for (i = i1; i <= i2; i = i + i3) {
                //
                if (work[i - 1] > vcrit) {
                    rec = one / vmax;
                    Rscal(n, rec, vr, 1);
                    Rscal(n, rec, vi, 1);
                    scale = scale * rec;
                    vmax = one;
                    vcrit = bignum;
                }
                //
                xr = vr[i - 1];
                xi = vi[i - 1];
                if (rightv) {
                    for (j = i + 1; j <= n; j = j + 1) {
                        xr = xr - b[(i - 1) + (j - 1) * ldb] * vr[j - 1] + b[((j + 1) - 1) + (i - 1) * ldb] * vi[j - 1];
                        xi = xi - b[(i - 1) + (j - 1) * ldb] * vi[j - 1] - b[((j + 1) - 1) + (i - 1) * ldb] * vr[j - 1];
                    }
                } else {
                    for (j = 1; j <= i - 1; j = j + 1) {
                        xr = xr - b[(j - 1) + (i - 1) * ldb] * vr[j - 1] + b[((i + 1) - 1) + (j - 1) * ldb] * vi[j - 1];
                        xi = xi - b[(j - 1) + (i - 1) * ldb] * vi[j - 1] - b[((i + 1) - 1) + (j - 1) * ldb] * vr[j - 1];
                    }
                }
                //
                w = abs(b[(i - 1) + (i - 1) * ldb]) + abs(b[((i + 1) - 1) + (i - 1) * ldb]);
                if (w > smlnum) {
                    if (w < one) {
                        w1 = abs(xr) + abs(xi);
                        if (w1 > w * bignum) {
                            rec = one / w1;
                            Rscal(n, rec, vr, 1);
                            Rscal(n, rec, vi, 1);
                            xr = vr[i - 1];
                            xi = vi[i - 1];
                            scale = scale * rec;
                            vmax = vmax * rec;
                        }
                    }
                    //
                    //                 Divide by diagonal element of B.
                    //
                    Rladiv(xr, xi, &b[(i - 1) + (i - 1) * ldb], &b[((i + 1) - 1) + (i - 1) * ldb], &vr[i - 1], vi[i - 1]);
                    vmax = max(abs(vr[i - 1]) + abs(vi[i - 1]), vmax);
                    vcrit = bignum / vmax;
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        vr[j - 1] = zero;
                        vi[j - 1] = zero;
                    }
                    vr[i - 1] = one;
                    vi[i - 1] = one;
                    scale = zero;
                    vmax = one;
                    vcrit = bignum;
                }
            }
            //
            //           Test for sufficient growth in the norm of (VR,VI).
            //
            vnorm = Rasum(n, vr, 1) + Rasum(n, vi, 1);
            if (vnorm >= growto * scale) {
                goto statement_280;
            }
            //
            //           Choose a new orthogonal starting vector and try again.
            //
            y = eps3 / (rootn + one);
            vr[1 - 1] = eps3;
            vi[1 - 1] = zero;
            //
            for (i = 2; i <= n; i = i + 1) {
                vr[i - 1] = y;
                vi[i - 1] = zero;
            }
            vr[(n - its + 1) - 1] = vr[(n - its + 1) - 1] - eps3 * rootn;
        }
        //
        //        Failure to find eigenvector in N iterations
        //
        info = 1;
    //
    statement_280:
        //
        //        Normalize eigenvector.
        //
        vnorm = zero;
        for (i = 1; i <= n; i = i + 1) {
            vnorm = max(vnorm, abs(vr[i - 1]) + abs(vi[i - 1]));
        }
        Rscal(n, one / vnorm, vr, 1);
        Rscal(n, one / vnorm, vi, 1);
        //
    }
    //
    //     End of Rlaein
    //
}

/*
 * Copyright (c) 2008-2022
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

void Ctgsy2(const char *trans, INTEGER const ijob, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *c, INTEGER const ldc, COMPLEX *d, INTEGER const ldd, COMPLEX *e, INTEGER const lde, COMPLEX *f, INTEGER const ldf, REAL &scale, REAL &rdsum, REAL &rdscal, INTEGER &info) {
    //
    //     Decode and test input parameters
    //
    info = 0;
    INTEGER ierr = 0;
    bool notran = Mlsame(trans, "N");
    if (!notran && !Mlsame(trans, "C")) {
        info = -1;
    } else if (notran) {
        if ((ijob < 0) || (ijob > 2)) {
            info = -2;
        }
    }
    if (info == 0) {
        if (m <= 0) {
            info = -3;
        } else if (n <= 0) {
            info = -4;
        } else if (lda < max((INTEGER)1, m)) {
            info = -6;
        } else if (ldb < max((INTEGER)1, n)) {
            info = -8;
        } else if (ldc < max((INTEGER)1, m)) {
            info = -10;
        } else if (ldd < max((INTEGER)1, m)) {
            info = -12;
        } else if (lde < max((INTEGER)1, n)) {
            info = -14;
        } else if (ldf < max((INTEGER)1, m)) {
            info = -16;
        }
    }
    if (info != 0) {
        Mxerbla("Ctgsy2", -info);
        return;
    }
    //
    const REAL one = 1.0;
    REAL scaloc = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    const INTEGER ldz = 2;
    COMPLEX z[ldz * ldz];
    COMPLEX rhs[ldz];
    INTEGER ipiv[ldz];
    INTEGER jpiv[ldz];
    INTEGER k = 0;
    const REAL zero = 0.0;
    COMPLEX alpha = 0.0;
    if (notran) {
        //
        //        Solve (I, J) - system
        //           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
        //           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
        //        for I = M, M - 1, ..., 1; J = 1, 2, ..., N
        //
        scale = one;
        scaloc = one;
        for (j = 1; j <= n; j = j + 1) {
            for (i = m; i >= 1; i = i - 1) {
                //
                //              Build 2 by 2 system
                //
                z[(1 - 1)] = a[(i - 1) + (i - 1) * lda];
                z[(2 - 1)] = d[(i - 1) + (i - 1) * ldd];
                z[(2 - 1) * ldz] = -b[(j - 1) + (j - 1) * ldb];
                z[(2 - 1) + (2 - 1) * ldz] = -e[(j - 1) + (j - 1) * lde];
                //
                //              Set up right hand side(s)
                //
                rhs[1 - 1] = c[(i - 1) + (j - 1) * ldc];
                rhs[2 - 1] = f[(i - 1) + (j - 1) * ldf];
                //
                //              Solve Z * x = RHS
                //
                Cgetc2(ldz, z, ldz, ipiv, jpiv, ierr);
                if (ierr > 0) {
                    info = ierr;
                }
                if (ijob == 0) {
                    Cgesc2(ldz, z, ldz, rhs, ipiv, jpiv, scaloc);
                    if (scaloc != one) {
                        for (k = 1; k <= n; k = k + 1) {
                            Cscal(m, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                            Cscal(m, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                        }
                        scale = scale * scaloc;
                    }
                } else {
                    Clatdf(ijob, ldz, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
                }
                //
                //              Unpack solution vector(s)
                //
                c[(i - 1) + (j - 1) * ldc] = rhs[1 - 1];
                f[(i - 1) + (j - 1) * ldf] = rhs[2 - 1];
                //
                //              Substitute R(I, J) and L(I, J) into remaining equation.
                //
                if (i > 1) {
                    alpha = -rhs[1 - 1];
                    Caxpy(i - 1, alpha, &a[(i - 1) * lda], 1, &c[(j - 1) * ldc], 1);
                    Caxpy(i - 1, alpha, &d[(i - 1) * ldd], 1, &f[(j - 1) * ldf], 1);
                }
                if (j < n) {
                    Caxpy(n - j, rhs[2 - 1], &b[(j - 1) + ((j + 1) - 1) * ldb], ldb, &c[(i - 1) + ((j + 1) - 1) * ldc], ldc);
                    Caxpy(n - j, rhs[2 - 1], &e[(j - 1) + ((j + 1) - 1) * lde], lde, &f[(i - 1) + ((j + 1) - 1) * ldf], ldf);
                }
                //
            }
        }
    } else {
        //
        //        Solve transposed (I, J) - system:
        //           A(I, I)**H * R(I, J) + D(I, I)**H * L(J, J) = C(I, J)
        //           R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
        //        for I = 1, 2, ..., M, J = N, N - 1, ..., 1
        //
        scale = one;
        scaloc = one;
        for (i = 1; i <= m; i = i + 1) {
            for (j = n; j >= 1; j = j - 1) {
                //
                //              Build 2 by 2 system Z**H
                //
                z[(1 - 1)] = conj(a[(i - 1) + (i - 1) * lda]);
                z[(2 - 1)] = -conj(b[(j - 1) + (j - 1) * ldb]);
                z[(2 - 1) * ldz] = conj(d[(i - 1) + (i - 1) * ldd]);
                z[(2 - 1) + (2 - 1) * ldz] = -conj(e[(j - 1) + (j - 1) * lde]);
                //
                //              Set up right hand side(s)
                //
                rhs[1 - 1] = c[(i - 1) + (j - 1) * ldc];
                rhs[2 - 1] = f[(i - 1) + (j - 1) * ldf];
                //
                //              Solve Z**H * x = RHS
                //
                Cgetc2(ldz, z, ldz, ipiv, jpiv, ierr);
                if (ierr > 0) {
                    info = ierr;
                }
                Cgesc2(ldz, z, ldz, rhs, ipiv, jpiv, scaloc);
                if (scaloc != one) {
                    for (k = 1; k <= n; k = k + 1) {
                        Cscal(m, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                        Cscal(m, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                    }
                    scale = scale * scaloc;
                }
                //
                //              Unpack solution vector(s)
                //
                c[(i - 1) + (j - 1) * ldc] = rhs[1 - 1];
                f[(i - 1) + (j - 1) * ldf] = rhs[2 - 1];
                //
                //              Substitute R(I, J) and L(I, J) into remaining equation.
                //
                for (k = 1; k <= j - 1; k = k + 1) {
                    f[(i - 1) + (k - 1) * ldf] += rhs[1 - 1] * conj(b[(k - 1) + (j - 1) * ldb]) + rhs[2 - 1] * conj(e[(k - 1) + (j - 1) * lde]);
                }
                for (k = i + 1; k <= m; k = k + 1) {
                    c[(k - 1) + (j - 1) * ldc] = c[(k - 1) + (j - 1) * ldc] - conj(a[(i - 1) + (k - 1) * lda]) * rhs[1 - 1] - conj(d[(i - 1) + (k - 1) * ldd]) * rhs[2 - 1];
                }
                //
            }
        }
    }
    //
    //     End of Ctgsy2
    //
}

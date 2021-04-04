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

void Cla_geamv(INTEGER const &trans, INTEGER const &m, INTEGER const &n, REAL const &alpha, COMPLEX *a, INTEGER const &lda, COMPLEX *x, INTEGER const &incx, REAL const &beta, REAL *y, INTEGER const &incy) {
    //
    //  -- LAPACK computational routine --
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function Definitions ..
    COMPLEX cdum = 0.0;
    abs1[cdum - 1] = abs(cdum.real()) + abs(cdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    if (!((trans == iMlatrans["N" - 1]) || (trans == iMlatrans["T" - 1]) || (trans == iMlatrans["C" - 1]))) {
        info = 1;
    } else if (m < 0) {
        info = 2;
    } else if (n < 0) {
        info = 3;
    } else if (lda < max((INTEGER)1, m)) {
        info = 6;
    } else if (incx == 0) {
        info = 8;
    } else if (incy == 0) {
        info = 11;
    }
    if (info != 0) {
        Mxerbla("Cla_geamv ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const COMPLEX zero = 0.0;
    const COMPLEX one = 1.0;
    if ((m == 0) || (n == 0) || ((alpha == zero) && (beta == one))) {
        return;
    }
    //
    //     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    //     up the start poINTEGERs in  X  and  Y.
    //
    INTEGER lenx = 0;
    INTEGER leny = 0;
    if (trans == iMlatrans["N" - 1]) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    INTEGER kx = 0;
    if (incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * incx;
    }
    INTEGER ky = 0;
    if (incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * incy;
    }
    //
    //     Set SAFE1 essentially to be the underflow threshold times the
    //     number of additions in each row.
    //
    REAL safe1 = dlamch("Safe minimum");
    safe1 = (n + 1) * safe1;
    //
    //     Form  y := alpha*abs(A)*abs(x) + beta*abs(y).
    //
    //     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
    //     the inexact flag.  Still doesn't help change the iteration order
    //     to per-column.
    //
    INTEGER iy = ky;
    INTEGER i = 0;
    bool symb_zero = false;
    INTEGER j = 0;
    REAL temp = 0.0;
    INTEGER jx = 0;
    if (incx == 1) {
        if (trans == iMlatrans["N" - 1]) {
            for (i = 1; i <= leny; i = i + 1) {
                if (beta == 0.0) {
                    symb_zero = true;
                    y[iy - 1] = 0.0;
                } else if (y[iy - 1] == 0.0) {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    y[iy - 1] = beta * abs(y[iy - 1]);
                }
                if (alpha != 0.0) {
                    for (j = 1; j <= lenx; j = j + 1) {
                        temp = abs1[a[(i - 1) + (j - 1) * lda] - 1];
                        symb_zero = symb_zero && (x[j - 1] == zero || temp == zero);
                        //
                        y[iy - 1] += alpha * abs1[x[j - 1] - 1] * temp;
                    }
                }
                //
                if (!symb_zero) {
                    y[iy - 1] += sign[(safe1 - 1) + (y[iy - 1] - 1) * ldsign];
                }
                //
                iy += incy;
            }
        } else {
            for (i = 1; i <= leny; i = i + 1) {
                if (beta == 0.0) {
                    symb_zero = true;
                    y[iy - 1] = 0.0;
                } else if (y[iy - 1] == 0.0) {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    y[iy - 1] = beta * abs(y[iy - 1]);
                }
                if (alpha != 0.0) {
                    for (j = 1; j <= lenx; j = j + 1) {
                        temp = abs1[a[(j - 1) + (i - 1) * lda] - 1];
                        symb_zero = symb_zero && (x[j - 1] == zero || temp == zero);
                        //
                        y[iy - 1] += alpha * abs1[x[j - 1] - 1] * temp;
                    }
                }
                //
                if (!symb_zero) {
                    y[iy - 1] += sign[(safe1 - 1) + (y[iy - 1] - 1) * ldsign];
                }
                //
                iy += incy;
            }
        }
    } else {
        if (trans == iMlatrans["N" - 1]) {
            for (i = 1; i <= leny; i = i + 1) {
                if (beta == 0.0) {
                    symb_zero = true;
                    y[iy - 1] = 0.0;
                } else if (y[iy - 1] == 0.0) {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    y[iy - 1] = beta * abs(y[iy - 1]);
                }
                if (alpha != 0.0) {
                    jx = kx;
                    for (j = 1; j <= lenx; j = j + 1) {
                        temp = abs1[a[(i - 1) + (j - 1) * lda] - 1];
                        symb_zero = symb_zero && (x[jx - 1] == zero || temp == zero);
                        //
                        y[iy - 1] += alpha * abs1[x[jx - 1] - 1] * temp;
                        jx += incx;
                    }
                }
                //
                if (!symb_zero) {
                    y[iy - 1] += sign[(safe1 - 1) + (y[iy - 1] - 1) * ldsign];
                }
                //
                iy += incy;
            }
        } else {
            for (i = 1; i <= leny; i = i + 1) {
                if (beta == 0.0) {
                    symb_zero = true;
                    y[iy - 1] = 0.0;
                } else if (y[iy - 1] == 0.0) {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    y[iy - 1] = beta * abs(y[iy - 1]);
                }
                if (alpha != 0.0) {
                    jx = kx;
                    for (j = 1; j <= lenx; j = j + 1) {
                        temp = abs1[a[(j - 1) + (i - 1) * lda] - 1];
                        symb_zero = symb_zero && (x[jx - 1] == zero || temp == zero);
                        //
                        y[iy - 1] += alpha * abs1[x[jx - 1] - 1] * temp;
                        jx += incx;
                    }
                }
                //
                if (!symb_zero) {
                    y[iy - 1] += sign[(safe1 - 1) + (y[iy - 1] - 1) * ldsign];
                }
                //
                iy += incy;
            }
        }
        //
    }
    //
    //     End of Cla_geamv
    //
}

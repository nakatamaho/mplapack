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

void Rlaqtr(bool const ltran, bool const lreal, INTEGER const n, REAL *t, INTEGER const ldt, REAL *b, REAL const w, REAL &scale, REAL *x, REAL *work, INTEGER &info) {
    bool notran = false;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL d[4];
    INTEGER ldd = 2;
    REAL xnorm = 0.0;
    REAL smin = 0.0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER n2 = 0;
    INTEGER n1 = 0;
    INTEGER k = 0;
    REAL xmax = 0.0;
    INTEGER jnext = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    REAL xj = 0.0;
    REAL tjj = 0.0;
    REAL tmp = 0.0;
    REAL rec = 0.0;
    REAL v[4];
    INTEGER ldv = 2;
    REAL scaloc = 0.0;
    INTEGER ierr = 0;
    REAL sminw = 0.0;
    REAL z = 0.0;
    REAL sr = 0.0;
    REAL si = 0.0;
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
    // =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Do not test the input parameters for errors
    //
    notran = !ltran;
    info = 0;
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Set constants to control overflow
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = one / smlnum;
    //
    xnorm = Rlange("M", n, n, t, ldt, d);
    if (!lreal) {
        xnorm = max({xnorm, abs(w), Rlange("M", n, 1, b, n, d)});
    }
    smin = max(smlnum, eps * xnorm);
    //
    //     Compute 1-norm of each column of strictly upper triangular
    //     part of T to control overflow in triangular solver.
    //
    work[1 - 1] = zero;
    for (j = 2; j <= n; j = j + 1) {
        work[j - 1] = Rasum(j - 1, &t[(j - 1) * ldt], 1);
    }
    //
    if (!lreal) {
        for (i = 2; i <= n; i = i + 1) {
            work[i - 1] += abs(b[i - 1]);
        }
    }
    //
    n2 = 2 * n;
    n1 = n;
    if (!lreal) {
        n1 = n2;
    }
    k = iRamax(n1, x, 1);
    xmax = abs(x[k - 1]);
    scale = one;
    //
    if (xmax > bignum) {
        scale = bignum / xmax;
        Rscal(n1, scale, x, 1);
        xmax = bignum;
    }
    //
    if (lreal) {
        //
        if (notran) {
            //
            //           Solve T*p = scale*c
            //
            jnext = n;
            for (j = n; j >= 1; j = j - 1) {
                if (j > jnext) {
                    goto statement_30;
                }
                j1 = j;
                j2 = j;
                jnext = j - 1;
                if (j > 1) {
                    if (t[(j - 1) + ((j - 1) - 1) * ldt] != zero) {
                        j1 = j - 1;
                        jnext = j - 2;
                    }
                }
                //
                if (j1 == j2) {
                    //
                    //                 Meet 1 by 1 diagonal block
                    //
                    //                 Scale to avoid overflow when computing
                    //                     x(j) = b(j)/T(j,j)
                    //
                    xj = abs(x[j1 - 1]);
                    tjj = abs(t[(j1 - 1) + (j1 - 1) * ldt]);
                    tmp = t[(j1 - 1) + (j1 - 1) * ldt];
                    if (tjj < smin) {
                        tmp = smin;
                        tjj = smin;
                        info = 1;
                    }
                    //
                    if (xj == zero) {
                        goto statement_30;
                    }
                    //
                    if (tjj < one) {
                        if (xj > bignum * tjj) {
                            rec = one / xj;
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    x[j1 - 1] = x[j1 - 1] / tmp;
                    xj = abs(x[j1 - 1]);
                    //
                    //                 Scale x if necessary to avoid overflow when adding a
                    //                 multiple of column j1 of T.
                    //
                    if (xj > one) {
                        rec = one / xj;
                        if (work[j1 - 1] > (bignum - xmax) * rec) {
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                        }
                    }
                    if (j1 > 1) {
                        Raxpy(j1 - 1, -x[j1 - 1], &t[(j1 - 1) * ldt], 1, x, 1);
                        k = iRamax(j1 - 1, x, 1);
                        xmax = abs(x[k - 1]);
                    }
                    //
                } else {
                    //
                    //                 Meet 2 by 2 diagonal block
                    //
                    //                 Call 2 by 2 linear system solve, to take
                    //                 care of possible overflow by scaling factor.
                    //
                    d[(1 - 1)] = x[j1 - 1];
                    d[(2 - 1)] = x[j2 - 1];
                    Rlaln2(false, 2, 1, smin, one, &t[(j1 - 1) + (j1 - 1) * ldt], ldt, one, one, d, 2, zero, zero, v, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 2;
                    }
                    //
                    if (scaloc != one) {
                        Rscal(n, scaloc, x, 1);
                        scale = scale * scaloc;
                    }
                    x[j1 - 1] = v[(1 - 1)];
                    x[j2 - 1] = v[(2 - 1)];
                    //
                    //                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
                    //                 to avoid overflow in updating right-hand side.
                    //
                    xj = max(abs(v[(1 - 1)]), abs(v[(2 - 1)]));
                    if (xj > one) {
                        rec = one / xj;
                        if (max(work[j1 - 1], work[j2 - 1]) > (bignum - xmax) * rec) {
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                        }
                    }
                    //
                    //                 Update right-hand side
                    //
                    if (j1 > 1) {
                        Raxpy(j1 - 1, -x[j1 - 1], &t[(j1 - 1) * ldt], 1, x, 1);
                        Raxpy(j1 - 1, -x[j2 - 1], &t[(j2 - 1) * ldt], 1, x, 1);
                        k = iRamax(j1 - 1, x, 1);
                        xmax = abs(x[k - 1]);
                    }
                    //
                }
            //
            statement_30:;
            }
            //
        } else {
            //
            //           Solve T**T*p = scale*c
            //
            jnext = 1;
            for (j = 1; j <= n; j = j + 1) {
                if (j < jnext) {
                    goto statement_40;
                }
                j1 = j;
                j2 = j;
                jnext = j + 1;
                if (j < n) {
                    if (t[((j + 1) - 1) + (j - 1) * ldt] != zero) {
                        j2 = j + 1;
                        jnext = j + 2;
                    }
                }
                //
                if (j1 == j2) {
                    //
                    //                 1 by 1 diagonal block
                    //
                    //                 Scale if necessary to avoid overflow in forming the
                    //                 right-hand side element by inner product.
                    //
                    xj = abs(x[j1 - 1]);
                    if (xmax > one) {
                        rec = one / xmax;
                        if (work[j1 - 1] > (bignum - xj) * rec) {
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    //
                    x[j1 - 1] = x[j1 - 1] - Rdot(j1 - 1, &t[(j1 - 1) * ldt], 1, x, 1);
                    //
                    xj = abs(x[j1 - 1]);
                    tjj = abs(t[(j1 - 1) + (j1 - 1) * ldt]);
                    tmp = t[(j1 - 1) + (j1 - 1) * ldt];
                    if (tjj < smin) {
                        tmp = smin;
                        tjj = smin;
                        info = 1;
                    }
                    //
                    if (tjj < one) {
                        if (xj > bignum * tjj) {
                            rec = one / xj;
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    x[j1 - 1] = x[j1 - 1] / tmp;
                    xmax = max(xmax, abs(x[j1 - 1]));
                    //
                } else {
                    //
                    //                 2 by 2 diagonal block
                    //
                    //                 Scale if necessary to avoid overflow in forming the
                    //                 right-hand side elements by inner product.
                    //
                    xj = max(abs(x[j1 - 1]), abs(x[j2 - 1]));
                    if (xmax > one) {
                        rec = one / xmax;
                        if (max(work[j2 - 1], work[j1 - 1]) > (bignum - xj) * rec) {
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    //
                    d[(1 - 1)] = x[j1 - 1] - Rdot(j1 - 1, &t[(j1 - 1) * ldt], 1, x, 1);
                    d[(2 - 1)] = x[j2 - 1] - Rdot(j1 - 1, &t[(j2 - 1) * ldt], 1, x, 1);
                    //
                    Rlaln2(true, 2, 1, smin, one, &t[(j1 - 1) + (j1 - 1) * ldt], ldt, one, one, d, 2, zero, zero, v, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 2;
                    }
                    //
                    if (scaloc != one) {
                        Rscal(n, scaloc, x, 1);
                        scale = scale * scaloc;
                    }
                    x[j1 - 1] = v[(1 - 1)];
                    x[j2 - 1] = v[(2 - 1)];
                    xmax = max({abs(x[j1 - 1]), abs(x[j2 - 1]), xmax});
                    //
                }
            statement_40:;
            }
        }
        //
    } else {
        //
        sminw = max(eps * abs(w), smin);
        if (notran) {
            //
            //           Solve (T + iB)*(p+iq) = c+id
            //
            jnext = n;
            for (j = n; j >= 1; j = j - 1) {
                if (j > jnext) {
                    goto statement_70;
                }
                j1 = j;
                j2 = j;
                jnext = j - 1;
                if (j > 1) {
                    if (t[(j - 1) + ((j - 1) - 1) * ldt] != zero) {
                        j1 = j - 1;
                        jnext = j - 2;
                    }
                }
                //
                if (j1 == j2) {
                    //
                    //                 1 by 1 diagonal block
                    //
                    //                 Scale if necessary to avoid overflow in division
                    //
                    z = w;
                    if (j1 == 1) {
                        z = b[1 - 1];
                    }
                    xj = abs(x[j1 - 1]) + abs(x[(n + j1) - 1]);
                    tjj = abs(t[(j1 - 1) + (j1 - 1) * ldt]) + abs(z);
                    tmp = t[(j1 - 1) + (j1 - 1) * ldt];
                    if (tjj < sminw) {
                        tmp = sminw;
                        tjj = sminw;
                        info = 1;
                    }
                    //
                    if (xj == zero) {
                        goto statement_70;
                    }
                    //
                    if (tjj < one) {
                        if (xj > bignum * tjj) {
                            rec = one / xj;
                            Rscal(n2, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    Rladiv(x[j1 - 1], x[(n + j1) - 1], tmp, z, sr, si);
                    x[j1 - 1] = sr;
                    x[(n + j1) - 1] = si;
                    xj = abs(x[j1 - 1]) + abs(x[(n + j1) - 1]);
                    //
                    //                 Scale x if necessary to avoid overflow when adding a
                    //                 multiple of column j1 of T.
                    //
                    if (xj > one) {
                        rec = one / xj;
                        if (work[j1 - 1] > (bignum - xmax) * rec) {
                            Rscal(n2, rec, x, 1);
                            scale = scale * rec;
                        }
                    }
                    //
                    if (j1 > 1) {
                        Raxpy(j1 - 1, -x[j1 - 1], &t[(j1 - 1) * ldt], 1, x, 1);
                        Raxpy(j1 - 1, -x[(n + j1) - 1], &t[(j1 - 1) * ldt], 1, &x[(n + 1) - 1], 1);
                        //
                        x[1 - 1] += b[j1 - 1] * x[(n + j1) - 1];
                        x[(n + 1) - 1] = x[(n + 1) - 1] - b[j1 - 1] * x[j1 - 1];
                        //
                        xmax = zero;
                        for (k = 1; k <= j1 - 1; k = k + 1) {
                            xmax = max(xmax, abs(x[k - 1]) + abs(x[(k + n) - 1]));
                        }
                    }
                    //
                } else {
                    //
                    //                 Meet 2 by 2 diagonal block
                    //
                    d[(1 - 1)] = x[j1 - 1];
                    d[(2 - 1)] = x[j2 - 1];
                    d[(2 - 1) * ldd] = x[(n + j1) - 1];
                    d[(2 - 1) + (2 - 1) * ldd] = x[(n + j2) - 1];
                    Rlaln2(false, 2, 2, sminw, one, &t[(j1 - 1) + (j1 - 1) * ldt], ldt, one, one, d, 2, zero, -w, v, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 2;
                    }
                    //
                    if (scaloc != one) {
                        Rscal(2 * n, scaloc, x, 1);
                        scale = scaloc * scale;
                    }
                    x[j1 - 1] = v[(1 - 1)];
                    x[j2 - 1] = v[(2 - 1)];
                    x[(n + j1) - 1] = v[(2 - 1) * ldv];
                    x[(n + j2) - 1] = v[(2 - 1) + (2 - 1) * ldv];
                    //
                    //                 Scale X(J1), .... to avoid overflow in
                    //                 updating right hand side.
                    //
                    xj = max(abs(v[(1 - 1)]) + abs(v[(2 - 1) * ldv]), abs(v[(2 - 1)]) + abs(v[(2 - 1) + (2 - 1) * ldv]));
                    if (xj > one) {
                        rec = one / xj;
                        if (max(work[j1 - 1], work[j2 - 1]) > (bignum - xmax) * rec) {
                            Rscal(n2, rec, x, 1);
                            scale = scale * rec;
                        }
                    }
                    //
                    //                 Update the right-hand side.
                    //
                    if (j1 > 1) {
                        Raxpy(j1 - 1, -x[j1 - 1], &t[(j1 - 1) * ldt], 1, x, 1);
                        Raxpy(j1 - 1, -x[j2 - 1], &t[(j2 - 1) * ldt], 1, x, 1);
                        //
                        Raxpy(j1 - 1, -x[(n + j1) - 1], &t[(j1 - 1) * ldt], 1, &x[(n + 1) - 1], 1);
                        Raxpy(j1 - 1, -x[(n + j2) - 1], &t[(j2 - 1) * ldt], 1, &x[(n + 1) - 1], 1);
                        //
                        x[1 - 1] += b[j1 - 1] * x[(n + j1) - 1] + b[j2 - 1] * x[(n + j2) - 1];
                        x[(n + 1) - 1] = x[(n + 1) - 1] - b[j1 - 1] * x[j1 - 1] - b[j2 - 1] * x[j2 - 1];
                        //
                        xmax = zero;
                        for (k = 1; k <= j1 - 1; k = k + 1) {
                            xmax = max(abs(x[k - 1]) + abs(x[(k + n) - 1]), xmax);
                        }
                    }
                    //
                }
            statement_70:;
            }
            //
        } else {
            //
            //           Solve (T + iB)**T*(p+iq) = c+id
            //
            jnext = 1;
            for (j = 1; j <= n; j = j + 1) {
                if (j < jnext) {
                    goto statement_80;
                }
                j1 = j;
                j2 = j;
                jnext = j + 1;
                if (j < n) {
                    if (t[((j + 1) - 1) + (j - 1) * ldt] != zero) {
                        j2 = j + 1;
                        jnext = j + 2;
                    }
                }
                //
                if (j1 == j2) {
                    //
                    //                 1 by 1 diagonal block
                    //
                    //                 Scale if necessary to avoid overflow in forming the
                    //                 right-hand side element by inner product.
                    //
                    xj = abs(x[j1 - 1]) + abs(x[(j1 + n) - 1]);
                    if (xmax > one) {
                        rec = one / xmax;
                        if (work[j1 - 1] > (bignum - xj) * rec) {
                            Rscal(n2, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    //
                    x[j1 - 1] = x[j1 - 1] - Rdot(j1 - 1, &t[(j1 - 1) * ldt], 1, x, 1);
                    x[(n + j1) - 1] = x[(n + j1) - 1] - Rdot(j1 - 1, &t[(j1 - 1) * ldt], 1, &x[(n + 1) - 1], 1);
                    if (j1 > 1) {
                        x[j1 - 1] = x[j1 - 1] - b[j1 - 1] * x[(n + 1) - 1];
                        x[(n + j1) - 1] += b[j1 - 1] * x[1 - 1];
                    }
                    xj = abs(x[j1 - 1]) + abs(x[(j1 + n) - 1]);
                    //
                    z = w;
                    if (j1 == 1) {
                        z = b[1 - 1];
                    }
                    //
                    //                 Scale if necessary to avoid overflow in
                    //                 complex division
                    //
                    tjj = abs(t[(j1 - 1) + (j1 - 1) * ldt]) + abs(z);
                    tmp = t[(j1 - 1) + (j1 - 1) * ldt];
                    if (tjj < sminw) {
                        tmp = sminw;
                        tjj = sminw;
                        info = 1;
                    }
                    //
                    if (tjj < one) {
                        if (xj > bignum * tjj) {
                            rec = one / xj;
                            Rscal(n2, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    Rladiv(x[j1 - 1], x[(n + j1) - 1], tmp, -z, sr, si);
                    x[j1 - 1] = sr;
                    x[(j1 + n) - 1] = si;
                    xmax = max(abs(x[j1 - 1]) + abs(x[(j1 + n) - 1]), xmax);
                    //
                } else {
                    //
                    //                 2 by 2 diagonal block
                    //
                    //                 Scale if necessary to avoid overflow in forming the
                    //                 right-hand side element by inner product.
                    //
                    xj = max(abs(x[j1 - 1]) + abs(x[(n + j1) - 1]), abs(x[j2 - 1]) + abs(x[(n + j2) - 1]));
                    if (xmax > one) {
                        rec = one / xmax;
                        if (max(work[j1 - 1], work[j2 - 1]) > (bignum - xj) / xmax) {
                            Rscal(n2, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    //
                    d[(1 - 1)] = x[j1 - 1] - Rdot(j1 - 1, &t[(j1 - 1) * ldt], 1, x, 1);
                    d[(2 - 1)] = x[j2 - 1] - Rdot(j1 - 1, &t[(j2 - 1) * ldt], 1, x, 1);
                    d[(2 - 1) * ldd] = x[(n + j1) - 1] - Rdot(j1 - 1, &t[(j1 - 1) * ldt], 1, &x[(n + 1) - 1], 1);
                    d[(2 - 1) + (2 - 1) * ldd] = x[(n + j2) - 1] - Rdot(j1 - 1, &t[(j2 - 1) * ldt], 1, &x[(n + 1) - 1], 1);
                    d[(1 - 1)] = d[(1 - 1)] - b[j1 - 1] * x[(n + 1) - 1];
                    d[(2 - 1)] = d[(2 - 1)] - b[j2 - 1] * x[(n + 1) - 1];
                    d[(2 - 1) * ldd] += b[j1 - 1] * x[1 - 1];
                    d[(2 - 1) + (2 - 1) * ldd] += b[j2 - 1] * x[1 - 1];
                    //
                    Rlaln2(true, 2, 2, sminw, one, &t[(j1 - 1) + (j1 - 1) * ldt], ldt, one, one, d, 2, zero, w, v, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 2;
                    }
                    //
                    if (scaloc != one) {
                        Rscal(n2, scaloc, x, 1);
                        scale = scaloc * scale;
                    }
                    x[j1 - 1] = v[(1 - 1)];
                    x[j2 - 1] = v[(2 - 1)];
                    x[(n + j1) - 1] = v[(2 - 1) * ldv];
                    x[(n + j2) - 1] = v[(2 - 1) + (2 - 1) * ldv];
                    xmax = max({abs(x[j1 - 1]) + abs(x[(n + j1) - 1]), abs(x[j2 - 1]) + abs(x[(n + j2) - 1]), xmax});
                    //
                }
            //
            statement_80:;
            }
            //
        }
        //
    }
    //
    //     End of Rlaqtr
    //
}

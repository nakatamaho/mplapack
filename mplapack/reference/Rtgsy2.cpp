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

void Rtgsy2(const char *trans, INTEGER const ijob, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *c, INTEGER const ldc, REAL *d, INTEGER const ldd, REAL *e, INTEGER const lde, REAL *f, INTEGER const ldf, REAL &scale, REAL const rdsum, REAL const rRscal, INTEGER *iwork, INTEGER &pq, INTEGER &info) {
    INTEGER ierr = 0;
    bool notran = false;
    INTEGER p = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    INTEGER q = 0;
    INTEGER j = 0;
    const REAL one = 1.0;
    REAL scaloc = 0.0;
    INTEGER js = 0;
    INTEGER jsp1 = 0;
    INTEGER je = 0;
    INTEGER nb = 0;
    INTEGER is = 0;
    INTEGER isp1 = 0;
    INTEGER ie = 0;
    INTEGER mb = 0;
    INTEGER zdim = 0;
    const INTEGER ldz = 8;
    arr_2d<ldz, ldz, REAL> z(fill0);
    arr_1d<ldz, REAL> rhs(fill0);
    arr_1d<ldz, int> ipiv(fill0);
    arr_1d<ldz, int> jpiv(fill0);
    INTEGER k = 0;
    REAL alpha = 0.0;
    INTEGER ii = 0;
    INTEGER jj = 0;
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
    //  Replaced various illegal calls to Rcopy by calls to Rlaset.
    //  Sven Hammarling, 27/5/02.
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
    //     Decode and test input parameters
    //
    info = 0;
    ierr = 0;
    notran = Mlsame(trans, "N");
    if (!notran && !Mlsame(trans, "T")) {
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
        Mxerbla("Rtgsy2", -info);
        return;
    }
    //
    //     Determine block structure of A
    //
    pq = 0;
    p = 0;
    i = 1;
statement_10:
    if (i > m) {
        goto statement_20;
    }
    p++;
    iwork[p - 1] = i;
    if (i == m) {
        goto statement_20;
    }
    if (a[((i + 1) - 1) + (i - 1) * lda] != zero) {
        i += 2;
    } else {
        i++;
    }
    goto statement_10;
statement_20:
    iwork[(p + 1) - 1] = m + 1;
    //
    //     Determine block structure of B
    //
    q = p + 1;
    j = 1;
statement_30:
    if (j > n) {
        goto statement_40;
    }
    q++;
    iwork[q - 1] = j;
    if (j == n) {
        goto statement_40;
    }
    if (b[((j + 1) - 1) + (j - 1) * ldb] != zero) {
        j += 2;
    } else {
        j++;
    }
    goto statement_30;
statement_40:
    iwork[(q + 1) - 1] = n + 1;
    pq = p * (q - p - 1);
    //
    if (notran) {
        //
        //        Solve (I, J) - subsystem
        //           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
        //           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
        //        for I = P, P - 1, ..., 1; J = 1, 2, ..., Q
        //
        scale = one;
        scaloc = one;
        for (j = p + 2; j <= q; j = j + 1) {
            js = iwork[j - 1];
            jsp1 = js + 1;
            je = iwork[(j + 1) - 1] - 1;
            nb = je - js + 1;
            for (i = p; i >= 1; i = i - 1) {
                //
                is = iwork[i - 1];
                isp1 = is + 1;
                ie = iwork[(i + 1) - 1] - 1;
                mb = ie - is + 1;
                zdim = mb * nb * 2;
                //
                if ((mb == 1) && (nb == 1)) {
                    //
                    //                 Build a 2-by-2 system Z * x = RHS
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = d[(is - 1) + (is - 1) * ldd];
                    z[(2 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    z[(2 - 1) + (2 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    rhs[1 - 1] = c[(is - 1) + (js - 1) * ldc];
                    rhs[2 - 1] = f[(is - 1) + (js - 1) * ldf];
                    //
                    //                 Solve Z * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    //
                    if (ijob == 0) {
                        Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                        if (scaloc != one) {
                            for (k = 1; k <= n; k = k + 1) {
                                Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                                Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                            }
                            scale = scale * scaloc;
                        }
                    } else {
                        Rlatdf(ijob, zdim, z, ldz, rhs, rdsum, rRscal, ipiv, jpiv);
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    c[(is - 1) + (js - 1) * ldc] = rhs[1 - 1];
                    f[(is - 1) + (js - 1) * ldf] = rhs[2 - 1];
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (i > 1) {
                        alpha = -rhs[1 - 1];
                        Raxpy(is - 1, alpha, &a[(is - 1) * lda], 1, &c[(js - 1) * ldc], 1);
                        Raxpy(is - 1, alpha, &d[(is - 1) * ldd], 1, f[(js - 1) * ldf], 1);
                    }
                    if (j < q) {
                        Raxpy(n - je, rhs[2 - 1], &b[(js - 1) + ((je + 1) - 1) * ldb], ldb, &c[(is - 1) + ((je + 1) - 1) * ldc], ldc);
                        Raxpy(n - je, rhs[2 - 1], &e[(js - 1) + ((je + 1) - 1) * lde], lde, f[(is - 1) + ((je + 1) - 1) * ldf], ldf);
                    }
                    //
                } else if ((mb == 1) && (nb == 2)) {
                    //
                    //                 Build a 4-by-4 system Z * x = RHS
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = zero;
                    z[(3 - 1)] = d[(is - 1) + (is - 1) * ldd];
                    z[(4 - 1)] = zero;
                    //
                    z[(2 - 1) * ldz] = zero;
                    z[(2 - 1) + (2 - 1) * ldz] = a[(is - 1) + (is - 1) * lda];
                    z[(3 - 1) + (2 - 1) * ldz] = zero;
                    z[(4 - 1) + (2 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    //
                    z[(3 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    z[(2 - 1) + (3 - 1) * ldz] = -b[(js - 1) + (jsp1 - 1) * ldb];
                    z[(3 - 1) + (3 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    z[(4 - 1) + (3 - 1) * ldz] = -e[(js - 1) + (jsp1 - 1) * lde];
                    //
                    z[(4 - 1) * ldz] = -b[(jsp1 - 1) + (js - 1) * ldb];
                    z[(2 - 1) + (4 - 1) * ldz] = -b[(jsp1 - 1) + (jsp1 - 1) * ldb];
                    z[(3 - 1) + (4 - 1) * ldz] = zero;
                    z[(4 - 1) + (4 - 1) * ldz] = -e[(jsp1 - 1) + (jsp1 - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    rhs[1 - 1] = c[(is - 1) + (js - 1) * ldc];
                    rhs[2 - 1] = c[(is - 1) + (jsp1 - 1) * ldc];
                    rhs[3 - 1] = f[(is - 1) + (js - 1) * ldf];
                    rhs[4 - 1] = f[(is - 1) + (jsp1 - 1) * ldf];
                    //
                    //                 Solve Z * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    //
                    if (ijob == 0) {
                        Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                        if (scaloc != one) {
                            for (k = 1; k <= n; k = k + 1) {
                                Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                                Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                            }
                            scale = scale * scaloc;
                        }
                    } else {
                        Rlatdf(ijob, zdim, z, ldz, rhs, rdsum, rRscal, ipiv, jpiv);
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    c[(is - 1) + (js - 1) * ldc] = rhs[1 - 1];
                    c[(is - 1) + (jsp1 - 1) * ldc] = rhs[2 - 1];
                    f[(is - 1) + (js - 1) * ldf] = rhs[3 - 1];
                    f[(is - 1) + (jsp1 - 1) * ldf] = rhs[4 - 1];
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (i > 1) {
                        Rger(is - 1, nb, -one, &a[(is - 1) * lda], 1, rhs[1 - 1], 1, &c[(js - 1) * ldc], ldc);
                        Rger(is - 1, nb, -one, &d[(is - 1) * ldd], 1, rhs[1 - 1], 1, f[(js - 1) * ldf], ldf);
                    }
                    if (j < q) {
                        Raxpy(n - je, rhs[3 - 1], &b[(js - 1) + ((je + 1) - 1) * ldb], ldb, &c[(is - 1) + ((je + 1) - 1) * ldc], ldc);
                        Raxpy(n - je, rhs[3 - 1], &e[(js - 1) + ((je + 1) - 1) * lde], lde, f[(is - 1) + ((je + 1) - 1) * ldf], ldf);
                        Raxpy(n - je, rhs[4 - 1], &b[(jsp1 - 1) + ((je + 1) - 1) * ldb], ldb, &c[(is - 1) + ((je + 1) - 1) * ldc], ldc);
                        Raxpy(n - je, rhs[4 - 1], &e[(jsp1 - 1) + ((je + 1) - 1) * lde], lde, f[(is - 1) + ((je + 1) - 1) * ldf], ldf);
                    }
                    //
                } else if ((mb == 2) && (nb == 1)) {
                    //
                    //                 Build a 4-by-4 system Z * x = RHS
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = a[(isp1 - 1) + (is - 1) * lda];
                    z[(3 - 1)] = d[(is - 1) + (is - 1) * ldd];
                    z[(4 - 1)] = zero;
                    //
                    z[(2 - 1) * ldz] = a[(is - 1) + (isp1 - 1) * lda];
                    z[(2 - 1) + (2 - 1) * ldz] = a[(isp1 - 1) + (isp1 - 1) * lda];
                    z[(3 - 1) + (2 - 1) * ldz] = d[(is - 1) + (isp1 - 1) * ldd];
                    z[(4 - 1) + (2 - 1) * ldz] = d[(isp1 - 1) + (isp1 - 1) * ldd];
                    //
                    z[(3 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    z[(2 - 1) + (3 - 1) * ldz] = zero;
                    z[(3 - 1) + (3 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    z[(4 - 1) + (3 - 1) * ldz] = zero;
                    //
                    z[(4 - 1) * ldz] = zero;
                    z[(2 - 1) + (4 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    z[(3 - 1) + (4 - 1) * ldz] = zero;
                    z[(4 - 1) + (4 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    rhs[1 - 1] = c[(is - 1) + (js - 1) * ldc];
                    rhs[2 - 1] = c[(isp1 - 1) + (js - 1) * ldc];
                    rhs[3 - 1] = f[(is - 1) + (js - 1) * ldf];
                    rhs[4 - 1] = f[(isp1 - 1) + (js - 1) * ldf];
                    //
                    //                 Solve Z * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    if (ijob == 0) {
                        Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                        if (scaloc != one) {
                            for (k = 1; k <= n; k = k + 1) {
                                Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                                Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                            }
                            scale = scale * scaloc;
                        }
                    } else {
                        Rlatdf(ijob, zdim, z, ldz, rhs, rdsum, rRscal, ipiv, jpiv);
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    c[(is - 1) + (js - 1) * ldc] = rhs[1 - 1];
                    c[(isp1 - 1) + (js - 1) * ldc] = rhs[2 - 1];
                    f[(is - 1) + (js - 1) * ldf] = rhs[3 - 1];
                    f[(isp1 - 1) + (js - 1) * ldf] = rhs[4 - 1];
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (i > 1) {
                        Rgemv("N", is - 1, mb, -one, &a[(is - 1) * lda], lda, rhs[1 - 1], 1, one, &c[(js - 1) * ldc], 1);
                        Rgemv("N", is - 1, mb, -one, &d[(is - 1) * ldd], ldd, rhs[1 - 1], 1, one, f[(js - 1) * ldf], 1);
                    }
                    if (j < q) {
                        Rger(mb, n - je, one, rhs[3 - 1], 1, &b[(js - 1) + ((je + 1) - 1) * ldb], ldb, &c[(is - 1) + ((je + 1) - 1) * ldc], ldc);
                        Rger(mb, n - je, one, rhs[3 - 1], 1, &e[(js - 1) + ((je + 1) - 1) * lde], lde, f[(is - 1) + ((je + 1) - 1) * ldf], ldf);
                    }
                    //
                } else if ((mb == 2) && (nb == 2)) {
                    //
                    //                 Build an 8-by-8 system Z * x = RHS
                    //
                    Rlaset("F", ldz, ldz, zero, zero, z, ldz);
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = a[(isp1 - 1) + (is - 1) * lda];
                    z[(5 - 1)] = d[(is - 1) + (is - 1) * ldd];
                    //
                    z[(2 - 1) * ldz] = a[(is - 1) + (isp1 - 1) * lda];
                    z[(2 - 1) + (2 - 1) * ldz] = a[(isp1 - 1) + (isp1 - 1) * lda];
                    z[(5 - 1) + (2 - 1) * ldz] = d[(is - 1) + (isp1 - 1) * ldd];
                    z[(6 - 1) + (2 - 1) * ldz] = d[(isp1 - 1) + (isp1 - 1) * ldd];
                    //
                    z[(3 - 1) + (3 - 1) * ldz] = a[(is - 1) + (is - 1) * lda];
                    z[(4 - 1) + (3 - 1) * ldz] = a[(isp1 - 1) + (is - 1) * lda];
                    z[(7 - 1) + (3 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    //
                    z[(3 - 1) + (4 - 1) * ldz] = a[(is - 1) + (isp1 - 1) * lda];
                    z[(4 - 1) + (4 - 1) * ldz] = a[(isp1 - 1) + (isp1 - 1) * lda];
                    z[(7 - 1) + (4 - 1) * ldz] = d[(is - 1) + (isp1 - 1) * ldd];
                    z[(8 - 1) + (4 - 1) * ldz] = d[(isp1 - 1) + (isp1 - 1) * ldd];
                    //
                    z[(5 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    z[(3 - 1) + (5 - 1) * ldz] = -b[(js - 1) + (jsp1 - 1) * ldb];
                    z[(5 - 1) + (5 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    z[(7 - 1) + (5 - 1) * ldz] = -e[(js - 1) + (jsp1 - 1) * lde];
                    //
                    z[(2 - 1) + (6 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    z[(4 - 1) + (6 - 1) * ldz] = -b[(js - 1) + (jsp1 - 1) * ldb];
                    z[(6 - 1) + (6 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    z[(8 - 1) + (6 - 1) * ldz] = -e[(js - 1) + (jsp1 - 1) * lde];
                    //
                    z[(7 - 1) * ldz] = -b[(jsp1 - 1) + (js - 1) * ldb];
                    z[(3 - 1) + (7 - 1) * ldz] = -b[(jsp1 - 1) + (jsp1 - 1) * ldb];
                    z[(7 - 1) + (7 - 1) * ldz] = -e[(jsp1 - 1) + (jsp1 - 1) * lde];
                    //
                    z[(2 - 1) + (8 - 1) * ldz] = -b[(jsp1 - 1) + (js - 1) * ldb];
                    z[(4 - 1) + (8 - 1) * ldz] = -b[(jsp1 - 1) + (jsp1 - 1) * ldb];
                    z[(8 - 1) + (8 - 1) * ldz] = -e[(jsp1 - 1) + (jsp1 - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    k = 1;
                    ii = mb * nb + 1;
                    for (jj = 0; jj <= nb - 1; jj = jj + 1) {
                        Rcopy(mb, &c[(is - 1) + ((js + jj) - 1) * ldc], 1, rhs[k - 1], 1);
                        Rcopy(mb, f[(is - 1) + ((js + jj) - 1) * ldf], 1, rhs[ii - 1], 1);
                        k += mb;
                        ii += mb;
                    }
                    //
                    //                 Solve Z * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    if (ijob == 0) {
                        Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                        if (scaloc != one) {
                            for (k = 1; k <= n; k = k + 1) {
                                Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                                Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                            }
                            scale = scale * scaloc;
                        }
                    } else {
                        Rlatdf(ijob, zdim, z, ldz, rhs, rdsum, rRscal, ipiv, jpiv);
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    k = 1;
                    ii = mb * nb + 1;
                    for (jj = 0; jj <= nb - 1; jj = jj + 1) {
                        Rcopy(mb, rhs[k - 1], 1, &c[(is - 1) + ((js + jj) - 1) * ldc], 1);
                        Rcopy(mb, rhs[ii - 1], 1, f[(is - 1) + ((js + jj) - 1) * ldf], 1);
                        k += mb;
                        ii += mb;
                    }
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (i > 1) {
                        Rgemm("N", "N", is - 1, nb, mb, -one, &a[(is - 1) * lda], lda, rhs[1 - 1], mb, one, &c[(js - 1) * ldc], ldc);
                        Rgemm("N", "N", is - 1, nb, mb, -one, &d[(is - 1) * ldd], ldd, rhs[1 - 1], mb, one, f[(js - 1) * ldf], ldf);
                    }
                    if (j < q) {
                        k = mb * nb + 1;
                        Rgemm("N", "N", mb, n - je, nb, one, rhs[k - 1], mb, &b[(js - 1) + ((je + 1) - 1) * ldb], ldb, one, &c[(is - 1) + ((je + 1) - 1) * ldc], ldc);
                        Rgemm("N", "N", mb, n - je, nb, one, rhs[k - 1], mb, &e[(js - 1) + ((je + 1) - 1) * lde], lde, one, f[(is - 1) + ((je + 1) - 1) * ldf], ldf);
                    }
                    //
                }
                //
            }
        }
    } else {
        //
        //        Solve (I, J) - subsystem
        //             A(I, I)**T * R(I, J) + D(I, I)**T * L(J, J)  =  C(I, J)
        //             R(I, I)  * B(J, J) + L(I, J)  * E(J, J)  = -F(I, J)
        //        for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1
        //
        scale = one;
        scaloc = one;
        for (i = 1; i <= p; i = i + 1) {
            //
            is = iwork[i - 1];
            isp1 = is + 1;
            ie = iwork[(i + 1) - 1] - 1;
            mb = ie - is + 1;
            for (j = q; j >= p + 2; j = j - 1) {
                //
                js = iwork[j - 1];
                jsp1 = js + 1;
                je = iwork[(j + 1) - 1] - 1;
                nb = je - js + 1;
                zdim = mb * nb * 2;
                if ((mb == 1) && (nb == 1)) {
                    //
                    //                 Build a 2-by-2 system Z**T * x = RHS
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = -b[(js - 1) + (js - 1) * ldb];
                    z[(2 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    z[(2 - 1) + (2 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    rhs[1 - 1] = c[(is - 1) + (js - 1) * ldc];
                    rhs[2 - 1] = f[(is - 1) + (js - 1) * ldf];
                    //
                    //                 Solve Z**T * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    //
                    Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                    if (scaloc != one) {
                        for (k = 1; k <= n; k = k + 1) {
                            Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                            Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                        }
                        scale = scale * scaloc;
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    c[(is - 1) + (js - 1) * ldc] = rhs[1 - 1];
                    f[(is - 1) + (js - 1) * ldf] = rhs[2 - 1];
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (j > p + 2) {
                        alpha = rhs[1 - 1];
                        Raxpy(js - 1, alpha, &b[(js - 1) * ldb], 1, f[(is - 1)], ldf);
                        alpha = rhs[2 - 1];
                        Raxpy(js - 1, alpha, &e[(js - 1) * lde], 1, f[(is - 1)], ldf);
                    }
                    if (i < p) {
                        alpha = -rhs[1 - 1];
                        Raxpy(m - ie, alpha, &a[(is - 1) + ((ie + 1) - 1) * lda], lda, &c[((ie + 1) - 1) + (js - 1) * ldc], 1);
                        alpha = -rhs[2 - 1];
                        Raxpy(m - ie, alpha, &d[(is - 1) + ((ie + 1) - 1) * ldd], ldd, &c[((ie + 1) - 1) + (js - 1) * ldc], 1);
                    }
                    //
                } else if ((mb == 1) && (nb == 2)) {
                    //
                    //                 Build a 4-by-4 system Z**T * x = RHS
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = zero;
                    z[(3 - 1)] = -b[(js - 1) + (js - 1) * ldb];
                    z[(4 - 1)] = -b[(jsp1 - 1) + (js - 1) * ldb];
                    //
                    z[(2 - 1) * ldz] = zero;
                    z[(2 - 1) + (2 - 1) * ldz] = a[(is - 1) + (is - 1) * lda];
                    z[(3 - 1) + (2 - 1) * ldz] = -b[(js - 1) + (jsp1 - 1) * ldb];
                    z[(4 - 1) + (2 - 1) * ldz] = -b[(jsp1 - 1) + (jsp1 - 1) * ldb];
                    //
                    z[(3 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    z[(2 - 1) + (3 - 1) * ldz] = zero;
                    z[(3 - 1) + (3 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    z[(4 - 1) + (3 - 1) * ldz] = zero;
                    //
                    z[(4 - 1) * ldz] = zero;
                    z[(2 - 1) + (4 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    z[(3 - 1) + (4 - 1) * ldz] = -e[(js - 1) + (jsp1 - 1) * lde];
                    z[(4 - 1) + (4 - 1) * ldz] = -e[(jsp1 - 1) + (jsp1 - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    rhs[1 - 1] = c[(is - 1) + (js - 1) * ldc];
                    rhs[2 - 1] = c[(is - 1) + (jsp1 - 1) * ldc];
                    rhs[3 - 1] = f[(is - 1) + (js - 1) * ldf];
                    rhs[4 - 1] = f[(is - 1) + (jsp1 - 1) * ldf];
                    //
                    //                 Solve Z**T * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                    if (scaloc != one) {
                        for (k = 1; k <= n; k = k + 1) {
                            Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                            Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                        }
                        scale = scale * scaloc;
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    c[(is - 1) + (js - 1) * ldc] = rhs[1 - 1];
                    c[(is - 1) + (jsp1 - 1) * ldc] = rhs[2 - 1];
                    f[(is - 1) + (js - 1) * ldf] = rhs[3 - 1];
                    f[(is - 1) + (jsp1 - 1) * ldf] = rhs[4 - 1];
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (j > p + 2) {
                        Raxpy(js - 1, rhs[1 - 1], &b[(js - 1) * ldb], 1, f[(is - 1)], ldf);
                        Raxpy(js - 1, rhs[2 - 1], &b[(jsp1 - 1) * ldb], 1, f[(is - 1)], ldf);
                        Raxpy(js - 1, rhs[3 - 1], &e[(js - 1) * lde], 1, f[(is - 1)], ldf);
                        Raxpy(js - 1, rhs[4 - 1], &e[(jsp1 - 1) * lde], 1, f[(is - 1)], ldf);
                    }
                    if (i < p) {
                        Rger(m - ie, nb, -one, &a[(is - 1) + ((ie + 1) - 1) * lda], lda, rhs[1 - 1], 1, &c[((ie + 1) - 1) + (js - 1) * ldc], ldc);
                        Rger(m - ie, nb, -one, &d[(is - 1) + ((ie + 1) - 1) * ldd], ldd, rhs[3 - 1], 1, &c[((ie + 1) - 1) + (js - 1) * ldc], ldc);
                    }
                    //
                } else if ((mb == 2) && (nb == 1)) {
                    //
                    //                 Build a 4-by-4 system Z**T * x = RHS
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = a[(is - 1) + (isp1 - 1) * lda];
                    z[(3 - 1)] = -b[(js - 1) + (js - 1) * ldb];
                    z[(4 - 1)] = zero;
                    //
                    z[(2 - 1) * ldz] = a[(isp1 - 1) + (is - 1) * lda];
                    z[(2 - 1) + (2 - 1) * ldz] = a[(isp1 - 1) + (isp1 - 1) * lda];
                    z[(3 - 1) + (2 - 1) * ldz] = zero;
                    z[(4 - 1) + (2 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    //
                    z[(3 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    z[(2 - 1) + (3 - 1) * ldz] = d[(is - 1) + (isp1 - 1) * ldd];
                    z[(3 - 1) + (3 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    z[(4 - 1) + (3 - 1) * ldz] = zero;
                    //
                    z[(4 - 1) * ldz] = zero;
                    z[(2 - 1) + (4 - 1) * ldz] = d[(isp1 - 1) + (isp1 - 1) * ldd];
                    z[(3 - 1) + (4 - 1) * ldz] = zero;
                    z[(4 - 1) + (4 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    rhs[1 - 1] = c[(is - 1) + (js - 1) * ldc];
                    rhs[2 - 1] = c[(isp1 - 1) + (js - 1) * ldc];
                    rhs[3 - 1] = f[(is - 1) + (js - 1) * ldf];
                    rhs[4 - 1] = f[(isp1 - 1) + (js - 1) * ldf];
                    //
                    //                 Solve Z**T * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    //
                    Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                    if (scaloc != one) {
                        for (k = 1; k <= n; k = k + 1) {
                            Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                            Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                        }
                        scale = scale * scaloc;
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    c[(is - 1) + (js - 1) * ldc] = rhs[1 - 1];
                    c[(isp1 - 1) + (js - 1) * ldc] = rhs[2 - 1];
                    f[(is - 1) + (js - 1) * ldf] = rhs[3 - 1];
                    f[(isp1 - 1) + (js - 1) * ldf] = rhs[4 - 1];
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (j > p + 2) {
                        Rger(mb, js - 1, one, rhs[1 - 1], 1, &b[(js - 1) * ldb], 1, f[(is - 1)], ldf);
                        Rger(mb, js - 1, one, rhs[3 - 1], 1, &e[(js - 1) * lde], 1, f[(is - 1)], ldf);
                    }
                    if (i < p) {
                        Rgemv("T", mb, m - ie, -one, &a[(is - 1) + ((ie + 1) - 1) * lda], lda, rhs[1 - 1], 1, one, &c[((ie + 1) - 1) + (js - 1) * ldc], 1);
                        Rgemv("T", mb, m - ie, -one, &d[(is - 1) + ((ie + 1) - 1) * ldd], ldd, rhs[3 - 1], 1, one, &c[((ie + 1) - 1) + (js - 1) * ldc], 1);
                    }
                    //
                } else if ((mb == 2) && (nb == 2)) {
                    //
                    //                 Build an 8-by-8 system Z**T * x = RHS
                    //
                    Rlaset("F", ldz, ldz, zero, zero, z, ldz);
                    //
                    z[(1 - 1)] = a[(is - 1) + (is - 1) * lda];
                    z[(2 - 1)] = a[(is - 1) + (isp1 - 1) * lda];
                    z[(5 - 1)] = -b[(js - 1) + (js - 1) * ldb];
                    z[(7 - 1)] = -b[(jsp1 - 1) + (js - 1) * ldb];
                    //
                    z[(2 - 1) * ldz] = a[(isp1 - 1) + (is - 1) * lda];
                    z[(2 - 1) + (2 - 1) * ldz] = a[(isp1 - 1) + (isp1 - 1) * lda];
                    z[(6 - 1) + (2 - 1) * ldz] = -b[(js - 1) + (js - 1) * ldb];
                    z[(8 - 1) + (2 - 1) * ldz] = -b[(jsp1 - 1) + (js - 1) * ldb];
                    //
                    z[(3 - 1) + (3 - 1) * ldz] = a[(is - 1) + (is - 1) * lda];
                    z[(4 - 1) + (3 - 1) * ldz] = a[(is - 1) + (isp1 - 1) * lda];
                    z[(5 - 1) + (3 - 1) * ldz] = -b[(js - 1) + (jsp1 - 1) * ldb];
                    z[(7 - 1) + (3 - 1) * ldz] = -b[(jsp1 - 1) + (jsp1 - 1) * ldb];
                    //
                    z[(3 - 1) + (4 - 1) * ldz] = a[(isp1 - 1) + (is - 1) * lda];
                    z[(4 - 1) + (4 - 1) * ldz] = a[(isp1 - 1) + (isp1 - 1) * lda];
                    z[(6 - 1) + (4 - 1) * ldz] = -b[(js - 1) + (jsp1 - 1) * ldb];
                    z[(8 - 1) + (4 - 1) * ldz] = -b[(jsp1 - 1) + (jsp1 - 1) * ldb];
                    //
                    z[(5 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    z[(2 - 1) + (5 - 1) * ldz] = d[(is - 1) + (isp1 - 1) * ldd];
                    z[(5 - 1) + (5 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    //
                    z[(2 - 1) + (6 - 1) * ldz] = d[(isp1 - 1) + (isp1 - 1) * ldd];
                    z[(6 - 1) + (6 - 1) * ldz] = -e[(js - 1) + (js - 1) * lde];
                    //
                    z[(3 - 1) + (7 - 1) * ldz] = d[(is - 1) + (is - 1) * ldd];
                    z[(4 - 1) + (7 - 1) * ldz] = d[(is - 1) + (isp1 - 1) * ldd];
                    z[(5 - 1) + (7 - 1) * ldz] = -e[(js - 1) + (jsp1 - 1) * lde];
                    z[(7 - 1) + (7 - 1) * ldz] = -e[(jsp1 - 1) + (jsp1 - 1) * lde];
                    //
                    z[(4 - 1) + (8 - 1) * ldz] = d[(isp1 - 1) + (isp1 - 1) * ldd];
                    z[(6 - 1) + (8 - 1) * ldz] = -e[(js - 1) + (jsp1 - 1) * lde];
                    z[(8 - 1) + (8 - 1) * ldz] = -e[(jsp1 - 1) + (jsp1 - 1) * lde];
                    //
                    //                 Set up right hand side(s)
                    //
                    k = 1;
                    ii = mb * nb + 1;
                    for (jj = 0; jj <= nb - 1; jj = jj + 1) {
                        Rcopy(mb, &c[(is - 1) + ((js + jj) - 1) * ldc], 1, rhs[k - 1], 1);
                        Rcopy(mb, f[(is - 1) + ((js + jj) - 1) * ldf], 1, rhs[ii - 1], 1);
                        k += mb;
                        ii += mb;
                    }
                    //
                    //                 Solve Z**T * x = RHS
                    //
                    Rgetc2(zdim, z, ldz, ipiv, jpiv, ierr);
                    if (ierr > 0) {
                        info = ierr;
                    }
                    //
                    Rgesc2(zdim, z, ldz, rhs, ipiv, jpiv, scaloc);
                    if (scaloc != one) {
                        for (k = 1; k <= n; k = k + 1) {
                            Rscal(m, scaloc, &c[(k - 1) * ldc], 1);
                            Rscal(m, scaloc, f[(k - 1) * ldf], 1);
                        }
                        scale = scale * scaloc;
                    }
                    //
                    //                 Unpack solution vector(s)
                    //
                    k = 1;
                    ii = mb * nb + 1;
                    for (jj = 0; jj <= nb - 1; jj = jj + 1) {
                        Rcopy(mb, rhs[k - 1], 1, &c[(is - 1) + ((js + jj) - 1) * ldc], 1);
                        Rcopy(mb, rhs[ii - 1], 1, f[(is - 1) + ((js + jj) - 1) * ldf], 1);
                        k += mb;
                        ii += mb;
                    }
                    //
                    //                 Substitute R(I, J) and L(I, J) into remaining
                    //                 equation.
                    //
                    if (j > p + 2) {
                        Rgemm("N", "T", mb, js - 1, nb, one, &c[(is - 1) + (js - 1) * ldc], ldc, &b[(js - 1) * ldb], ldb, one, f[(is - 1)], ldf);
                        Rgemm("N", "T", mb, js - 1, nb, one, f[(is - 1) + (js - 1) * ldf], ldf, &e[(js - 1) * lde], lde, one, f[(is - 1)], ldf);
                    }
                    if (i < p) {
                        Rgemm("T", "N", m - ie, nb, mb, -one, &a[(is - 1) + ((ie + 1) - 1) * lda], lda, &c[(is - 1) + (js - 1) * ldc], ldc, one, &c[((ie + 1) - 1) + (js - 1) * ldc], ldc);
                        Rgemm("T", "N", m - ie, nb, mb, -one, &d[(is - 1) + ((ie + 1) - 1) * ldd], ldd, f[(is - 1) + (js - 1) * ldf], ldf, one, &c[((ie + 1) - 1) + (js - 1) * ldc], ldc);
                    }
                    //
                }
                //
            }
        }
        //
    }
    //
    //     End of Rtgsy2
    //
}

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

void Cunbdb(const char *trans, const char *signs, INTEGER const m, INTEGER const p, INTEGER const q, COMPLEX *x11, INTEGER const ldx11, COMPLEX *x12, INTEGER const ldx12, COMPLEX *x21, INTEGER const ldx21, COMPLEX *x22, INTEGER const ldx22, REAL *theta, REAL *phi, COMPLEX *taup1, COMPLEX *taup2, COMPLEX *tauq1, COMPLEX *tauq2, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    //  ====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions
    //     ..
    //     .. Executable Statements ..
    //
    //     Test input arguments
    //
    info = 0;
    bool colmajor = !Mlsame(trans, "T");
    const REAL realone = 1.0;
    REAL z1 = 0.0;
    REAL z2 = 0.0;
    REAL z3 = 0.0;
    REAL z4 = 0.0;
    if (!Mlsame(signs, "O")) {
        z1 = realone;
        z2 = realone;
        z3 = realone;
        z4 = realone;
    } else {
        z1 = realone;
        z2 = -realone;
        z3 = realone;
        z4 = -realone;
    }
    bool lquery = lwork == -1;
    //
    if (m < 0) {
        info = -3;
    } else if (p < 0 || p > m) {
        info = -4;
    } else if (q < 0 || q > p || q > m - p || q > m - q) {
        info = -5;
    } else if (colmajor && ldx11 < max((INTEGER)1, p)) {
        info = -7;
    } else if (!colmajor && ldx11 < max((INTEGER)1, q)) {
        info = -7;
    } else if (colmajor && ldx12 < max((INTEGER)1, p)) {
        info = -9;
    } else if (!colmajor && ldx12 < max((INTEGER)1, m - q)) {
        info = -9;
    } else if (colmajor && ldx21 < max((INTEGER)1, m - p)) {
        info = -11;
    } else if (!colmajor && ldx21 < max((INTEGER)1, q)) {
        info = -11;
    } else if (colmajor && ldx22 < max((INTEGER)1, m - p)) {
        info = -13;
    } else if (!colmajor && ldx22 < max((INTEGER)1, m - q)) {
        info = -13;
    }
    //
    //     Compute workspace
    //
    INTEGER lworkopt = 0;
    INTEGER lworkmin = 0;
    if (info == 0) {
        lworkopt = m - q;
        lworkmin = m - q;
        work[1 - 1] = lworkopt;
        if (lwork < lworkmin && !lquery) {
            info = -21;
        }
    }
    if (info != 0) {
        Mxerbla("xORBDB", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Handle column-major and row-major separately
    //
    INTEGER i = 0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    if (colmajor) {
        //
        //        Reduce columns 1, ..., Q of X11, X12, X21, and X22
        //
        for (i = 1; i <= q; i = i + 1) {
            //
            if (i == 1) {
                Cscal(p - i + 1, COMPLEX(z1, 0.0), &x11[(i - 1) + (i - 1) * ldx11], 1);
            } else {
                Cscal(p - i + 1, COMPLEX(z1 * cos(phi[(i - 1) - 1]), 0.0), &x11[(i - 1) + (i - 1) * ldx11], 1);
                Caxpy(p - i + 1, COMPLEX(-z1 * z3 * z4 * sin(phi[(i - 1) - 1]), 0.0), &x12[(i - 1) + ((i - 1) - 1) * ldx12], 1, &x11[(i - 1) + (i - 1) * ldx11], 1);
            }
            if (i == 1) {
                Cscal(m - p - i + 1, COMPLEX(z2, 0.0), &x21[(i - 1) + (i - 1) * ldx21], 1);
            } else {
                Cscal(m - p - i + 1, COMPLEX(z2 * cos(phi[(i - 1) - 1]), 0.0), &x21[(i - 1) + (i - 1) * ldx21], 1);
                Caxpy(m - p - i + 1, COMPLEX(-z2 * z3 * z4 * sin(phi[(i - 1) - 1]), 0.0), &x22[(i - 1) + ((i - 1) - 1) * ldx22], 1, &x21[(i - 1) + (i - 1) * ldx21], 1);
            }
            //
            theta[i - 1] = atan2(RCnrm2(m - p - i + 1, &x21[(i - 1) + (i - 1) * ldx21], 1), RCnrm2(p - i + 1, &x11[(i - 1) + (i - 1) * ldx11], 1));
            //
            if (p > i) {
                Clarfgp(p - i + 1, x11[(i - 1) + (i - 1) * ldx11], &x11[((i + 1) - 1) + (i - 1) * ldx11], 1, taup1[i - 1]);
            } else if (p == i) {
                Clarfgp(p - i + 1, x11[(i - 1) + (i - 1) * ldx11], &x11[(i - 1) + (i - 1) * ldx11], 1, taup1[i - 1]);
            }
            x11[(i - 1) + (i - 1) * ldx11] = one;
            if (m - p > i) {
                Clarfgp(m - p - i + 1, x21[(i - 1) + (i - 1) * ldx21], &x21[((i + 1) - 1) + (i - 1) * ldx21], 1, taup2[i - 1]);
            } else if (m - p == i) {
                Clarfgp(m - p - i + 1, x21[(i - 1) + (i - 1) * ldx21], &x21[(i - 1) + (i - 1) * ldx21], 1, taup2[i - 1]);
            }
            x21[(i - 1) + (i - 1) * ldx21] = one;
            //
            if (q > i) {
                Clarf("L", p - i + 1, q - i, &x11[(i - 1) + (i - 1) * ldx11], 1, conj(taup1[i - 1]), &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11, work);
                Clarf("L", m - p - i + 1, q - i, &x21[(i - 1) + (i - 1) * ldx21], 1, conj(taup2[i - 1]), &x21[(i - 1) + ((i + 1) - 1) * ldx21], ldx21, work);
            }
            if (m - q + 1 > i) {
                Clarf("L", p - i + 1, m - q - i + 1, &x11[(i - 1) + (i - 1) * ldx11], 1, conj(taup1[i - 1]), &x12[(i - 1) + (i - 1) * ldx12], ldx12, work);
                Clarf("L", m - p - i + 1, m - q - i + 1, &x21[(i - 1) + (i - 1) * ldx21], 1, conj(taup2[i - 1]), &x22[(i - 1) + (i - 1) * ldx22], ldx22, work);
            }
            //
            if (i < q) {
                Cscal(q - i, COMPLEX(-z1 * z3 * sin(theta[i - 1]), 0.0), &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11);
                Caxpy(q - i, COMPLEX(z2 * z3 * cos(theta[i - 1]), 0.0), &x21[(i - 1) + ((i + 1) - 1) * ldx21], ldx21, &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11);
            }
            Cscal(m - q - i + 1, COMPLEX(-z1 * z4 * sin(theta[i - 1]), 0.0), &x12[(i - 1) + (i - 1) * ldx12], ldx12);
            Caxpy(m - q - i + 1, COMPLEX(z2 * z4 * cos(theta[i - 1]), 0.0), &x22[(i - 1) + (i - 1) * ldx22], ldx22, &x12[(i - 1) + (i - 1) * ldx12], ldx12);
            //
            if (i < q) {
                phi[i - 1] = atan2(RCnrm2(q - i, &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11), RCnrm2(m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12));
            }
            //
            if (i < q) {
                Clacgv(q - i, &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11);
                if (i == q - 1) {
                    Clarfgp(q - i, x11[(i - 1) + ((i + 1) - 1) * ldx11], &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11, tauq1[i - 1]);
                } else {
                    Clarfgp(q - i, x11[(i - 1) + ((i + 1) - 1) * ldx11], &x11[(i - 1) + ((i + 2) - 1) * ldx11], ldx11, tauq1[i - 1]);
                }
                x11[(i - 1) + ((i + 1) - 1) * ldx11] = one;
            }
            if (m - q + 1 > i) {
                Clacgv(m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12);
                if (m - q == i) {
                    Clarfgp(m - q - i + 1, x12[(i - 1) + (i - 1) * ldx12], &x12[(i - 1) + (i - 1) * ldx12], ldx12, tauq2[i - 1]);
                } else {
                    Clarfgp(m - q - i + 1, x12[(i - 1) + (i - 1) * ldx12], &x12[(i - 1) + ((i + 1) - 1) * ldx12], ldx12, tauq2[i - 1]);
                }
            }
            x12[(i - 1) + (i - 1) * ldx12] = one;
            //
            if (i < q) {
                Clarf("R", p - i, q - i, &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11, tauq1[i - 1], &x11[((i + 1) - 1) + ((i + 1) - 1) * ldx11], ldx11, work);
                Clarf("R", m - p - i, q - i, &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11, tauq1[i - 1], &x21[((i + 1) - 1) + ((i + 1) - 1) * ldx21], ldx21, work);
            }
            if (p > i) {
                Clarf("R", p - i, m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12, tauq2[i - 1], &x12[((i + 1) - 1) + (i - 1) * ldx12], ldx12, work);
            }
            if (m - p > i) {
                Clarf("R", m - p - i, m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12, tauq2[i - 1], &x22[((i + 1) - 1) + (i - 1) * ldx22], ldx22, work);
            }
            //
            if (i < q) {
                Clacgv(q - i, &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11);
            }
            Clacgv(m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12);
            //
        }
        //
        //        Reduce columns Q + 1, ..., P of X12, X22
        //
        for (i = q + 1; i <= p; i = i + 1) {
            //
            Cscal(m - q - i + 1, COMPLEX(-z1 * z4, 0.0), &x12[(i - 1) + (i - 1) * ldx12], ldx12);
            Clacgv(m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12);
            if (i >= m - q) {
                Clarfgp(m - q - i + 1, x12[(i - 1) + (i - 1) * ldx12], &x12[(i - 1) + (i - 1) * ldx12], ldx12, tauq2[i - 1]);
            } else {
                Clarfgp(m - q - i + 1, x12[(i - 1) + (i - 1) * ldx12], &x12[(i - 1) + ((i + 1) - 1) * ldx12], ldx12, tauq2[i - 1]);
            }
            x12[(i - 1) + (i - 1) * ldx12] = one;
            //
            if (p > i) {
                Clarf("R", p - i, m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12, tauq2[i - 1], &x12[((i + 1) - 1) + (i - 1) * ldx12], ldx12, work);
            }
            if (m - p - q >= 1) {
                Clarf("R", m - p - q, m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12, tauq2[i - 1], &x22[((q + 1) - 1) + (i - 1) * ldx22], ldx22, work);
            }
            //
            Clacgv(m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], ldx12);
            //
        }
        //
        //        Reduce columns P + 1, ..., M - Q of X12, X22
        //
        for (i = 1; i <= m - p - q; i = i + 1) {
            //
            Cscal(m - p - q - i + 1, COMPLEX(z2 * z4, 0.0), &x22[((q + i) - 1) + ((p + i) - 1) * ldx22], ldx22);
            Clacgv(m - p - q - i + 1, &x22[((q + i) - 1) + ((p + i) - 1) * ldx22], ldx22);
            Clarfgp(m - p - q - i + 1, x22[((q + i) - 1) + ((p + i) - 1) * ldx22], &x22[((q + i) - 1) + ((p + i + 1) - 1) * ldx22], ldx22, tauq2[(p + i) - 1]);
            x22[((q + i) - 1) + ((p + i) - 1) * ldx22] = one;
            Clarf("R", m - p - q - i, m - p - q - i + 1, &x22[((q + i) - 1) + ((p + i) - 1) * ldx22], ldx22, tauq2[(p + i) - 1], &x22[((q + i + 1) - 1) + ((p + i) - 1) * ldx22], ldx22, work);
            //
            Clacgv(m - p - q - i + 1, &x22[((q + i) - 1) + ((p + i) - 1) * ldx22], ldx22);
            //
        }
        //
    } else {
        //
        //        Reduce columns 1, ..., Q of X11, X12, X21, X22
        //
        for (i = 1; i <= q; i = i + 1) {
            //
            if (i == 1) {
                Cscal(p - i + 1, COMPLEX(z1, 0.0), &x11[(i - 1) + (i - 1) * ldx11], ldx11);
            } else {
                Cscal(p - i + 1, COMPLEX(z1 * cos(phi[(i - 1) - 1]), 0.0), &x11[(i - 1) + (i - 1) * ldx11], ldx11);
                Caxpy(p - i + 1, COMPLEX(-z1 * z3 * z4 * sin(phi[(i - 1) - 1]), 0.0), &x12[((i - 1) - 1) + (i - 1) * ldx12], ldx12, &x11[(i - 1) + (i - 1) * ldx11], ldx11);
            }
            if (i == 1) {
                Cscal(m - p - i + 1, COMPLEX(z2, 0.0), &x21[(i - 1) + (i - 1) * ldx21], ldx21);
            } else {
                Cscal(m - p - i + 1, COMPLEX(z2 * cos(phi[(i - 1) - 1]), 0.0), &x21[(i - 1) + (i - 1) * ldx21], ldx21);
                Caxpy(m - p - i + 1, COMPLEX(-z2 * z3 * z4 * sin(phi[(i - 1) - 1]), 0.0), &x22[((i - 1) - 1) + (i - 1) * ldx22], ldx22, &x21[(i - 1) + (i - 1) * ldx21], ldx21);
            }
            //
            theta[i - 1] = atan2(RCnrm2(m - p - i + 1, &x21[(i - 1) + (i - 1) * ldx21], ldx21), RCnrm2(p - i + 1, &x11[(i - 1) + (i - 1) * ldx11], ldx11));
            //
            Clacgv(p - i + 1, &x11[(i - 1) + (i - 1) * ldx11], ldx11);
            Clacgv(m - p - i + 1, &x21[(i - 1) + (i - 1) * ldx21], ldx21);
            //
            Clarfgp(p - i + 1, x11[(i - 1) + (i - 1) * ldx11], &x11[(i - 1) + ((i + 1) - 1) * ldx11], ldx11, taup1[i - 1]);
            x11[(i - 1) + (i - 1) * ldx11] = one;
            if (i == m - p) {
                Clarfgp(m - p - i + 1, x21[(i - 1) + (i - 1) * ldx21], &x21[(i - 1) + (i - 1) * ldx21], ldx21, taup2[i - 1]);
            } else {
                Clarfgp(m - p - i + 1, x21[(i - 1) + (i - 1) * ldx21], &x21[(i - 1) + ((i + 1) - 1) * ldx21], ldx21, taup2[i - 1]);
            }
            x21[(i - 1) + (i - 1) * ldx21] = one;
            //
            Clarf("R", q - i, p - i + 1, &x11[(i - 1) + (i - 1) * ldx11], ldx11, taup1[i - 1], &x11[((i + 1) - 1) + (i - 1) * ldx11], ldx11, work);
            Clarf("R", m - q - i + 1, p - i + 1, &x11[(i - 1) + (i - 1) * ldx11], ldx11, taup1[i - 1], &x12[(i - 1) + (i - 1) * ldx12], ldx12, work);
            Clarf("R", q - i, m - p - i + 1, &x21[(i - 1) + (i - 1) * ldx21], ldx21, taup2[i - 1], &x21[((i + 1) - 1) + (i - 1) * ldx21], ldx21, work);
            Clarf("R", m - q - i + 1, m - p - i + 1, &x21[(i - 1) + (i - 1) * ldx21], ldx21, taup2[i - 1], &x22[(i - 1) + (i - 1) * ldx22], ldx22, work);
            //
            Clacgv(p - i + 1, &x11[(i - 1) + (i - 1) * ldx11], ldx11);
            Clacgv(m - p - i + 1, &x21[(i - 1) + (i - 1) * ldx21], ldx21);
            //
            if (i < q) {
                Cscal(q - i, COMPLEX(-z1 * z3 * sin(theta[i - 1]), 0.0), &x11[((i + 1) - 1) + (i - 1) * ldx11], 1);
                Caxpy(q - i, COMPLEX(z2 * z3 * cos(theta[i - 1]), 0.0), &x21[((i + 1) - 1) + (i - 1) * ldx21], 1, &x11[((i + 1) - 1) + (i - 1) * ldx11], 1);
            }
            Cscal(m - q - i + 1, COMPLEX(-z1 * z4 * sin(theta[i - 1]), 0.0), &x12[(i - 1) + (i - 1) * ldx12], 1);
            Caxpy(m - q - i + 1, COMPLEX(z2 * z4 * cos(theta[i - 1]), 0.0), &x22[(i - 1) + (i - 1) * ldx22], 1, &x12[(i - 1) + (i - 1) * ldx12], 1);
            //
            if (i < q) {
                phi[i - 1] = atan2(RCnrm2(q - i, &x11[((i + 1) - 1) + (i - 1) * ldx11], 1), RCnrm2(m - q - i + 1, &x12[(i - 1) + (i - 1) * ldx12], 1));
            }
            //
            if (i < q) {
                Clarfgp(q - i, x11[((i + 1) - 1) + (i - 1) * ldx11], &x11[((i + 2) - 1) + (i - 1) * ldx11], 1, tauq1[i - 1]);
                x11[((i + 1) - 1) + (i - 1) * ldx11] = one;
            }
            Clarfgp(m - q - i + 1, x12[(i - 1) + (i - 1) * ldx12], &x12[((i + 1) - 1) + (i - 1) * ldx12], 1, tauq2[i - 1]);
            x12[(i - 1) + (i - 1) * ldx12] = one;
            //
            if (i < q) {
                Clarf("L", q - i, p - i, &x11[((i + 1) - 1) + (i - 1) * ldx11], 1, conj(tauq1[i - 1]), &x11[((i + 1) - 1) + ((i + 1) - 1) * ldx11], ldx11, work);
                Clarf("L", q - i, m - p - i, &x11[((i + 1) - 1) + (i - 1) * ldx11], 1, conj(tauq1[i - 1]), &x21[((i + 1) - 1) + ((i + 1) - 1) * ldx21], ldx21, work);
            }
            Clarf("L", m - q - i + 1, p - i, &x12[(i - 1) + (i - 1) * ldx12], 1, conj(tauq2[i - 1]), &x12[(i - 1) + ((i + 1) - 1) * ldx12], ldx12, work);
            if (m - p > i) {
                Clarf("L", m - q - i + 1, m - p - i, &x12[(i - 1) + (i - 1) * ldx12], 1, conj(tauq2[i - 1]), &x22[(i - 1) + ((i + 1) - 1) * ldx22], ldx22, work);
            }
            //
        }
        //
        //        Reduce columns Q + 1, ..., P of X12, X22
        //
        for (i = q + 1; i <= p; i = i + 1) {
            //
            Cscal(m - q - i + 1, COMPLEX(-z1 * z4, 0.0), &x12[(i - 1) + (i - 1) * ldx12], 1);
            Clarfgp(m - q - i + 1, x12[(i - 1) + (i - 1) * ldx12], &x12[((i + 1) - 1) + (i - 1) * ldx12], 1, tauq2[i - 1]);
            x12[(i - 1) + (i - 1) * ldx12] = one;
            //
            if (p > i) {
                Clarf("L", m - q - i + 1, p - i, &x12[(i - 1) + (i - 1) * ldx12], 1, conj(tauq2[i - 1]), &x12[(i - 1) + ((i + 1) - 1) * ldx12], ldx12, work);
            }
            if (m - p - q >= 1) {
                Clarf("L", m - q - i + 1, m - p - q, &x12[(i - 1) + (i - 1) * ldx12], 1, conj(tauq2[i - 1]), &x22[(i - 1) + ((q + 1) - 1) * ldx22], ldx22, work);
            }
            //
        }
        //
        //        Reduce columns P + 1, ..., M - Q of X12, X22
        //
        for (i = 1; i <= m - p - q; i = i + 1) {
            //
            Cscal(m - p - q - i + 1, COMPLEX(z2 * z4, 0.0), &x22[((p + i) - 1) + ((q + i) - 1) * ldx22], 1);
            Clarfgp(m - p - q - i + 1, x22[((p + i) - 1) + ((q + i) - 1) * ldx22], &x22[((p + i + 1) - 1) + ((q + i) - 1) * ldx22], 1, tauq2[(p + i) - 1]);
            x22[((p + i) - 1) + ((q + i) - 1) * ldx22] = one;
            //
            if (m - p - q != i) {
                Clarf("L", m - p - q - i + 1, m - p - q - i, &x22[((p + i) - 1) + ((q + i) - 1) * ldx22], 1, conj(tauq2[(p + i) - 1]), &x22[((p + i) - 1) + ((q + i + 1) - 1) * ldx22], ldx22, work);
            }
            //
        }
        //
    }
    //
    //     End of Cunbdb
    //
}

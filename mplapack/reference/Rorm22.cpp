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

void Rorm22(const char *side, const char *trans, INTEGER const m, INTEGER const n, INTEGER const n1, INTEGER const n2, REAL *q, INTEGER const ldq, REAL *c, INTEGER const ldc, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //
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
    //     Test the input arguments
    //
    info = 0;
    bool left = Mlsame(side, "L");
    bool notran = Mlsame(trans, "N");
    bool lquery = (lwork == -1);
    //
    //     NQ is the order of Q;
    //
    INTEGER nq = 0;
    if (left) {
        nq = m;
    } else {
        nq = n;
    }
    INTEGER nw = nq;
    if (n1 == 0 || n2 == 0) {
        nw = 1;
    }
    if (!left && !Mlsame(side, "R")) {
        info = -1;
    } else if (!Mlsame(trans, "N") && !Mlsame(trans, "T")) {
        info = -2;
    } else if (m < 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (n1 < 0 || n1 + n2 != nq) {
        info = -5;
    } else if (n2 < 0) {
        info = -6;
    } else if (ldq < max((INTEGER)1, nq)) {
        info = -8;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -10;
    } else if (lwork < nw && !lquery) {
        info = -12;
    }
    //
    INTEGER lwkopt = 0;
    if (info == 0) {
        lwkopt = m * n;
        work[1 - 1] = castREAL(lwkopt);
    }
    //
    if (info != 0) {
        Mxerbla("Rorm22", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    //     Degenerate cases (N1 = 0 or N2 = 0) are handled using Rtrmm.
    //
    const REAL one = 1.0;
    if (n1 == 0) {
        Rtrmm(side, "Upper", trans, "Non-Unit", m, n, one, q, ldq, c, ldc);
        work[1 - 1] = one;
        return;
    } else if (n2 == 0) {
        Rtrmm(side, "Lower", trans, "Non-Unit", m, n, one, q, ldq, c, ldc);
        work[1 - 1] = one;
        return;
    }
    //
    //     Compute the largest chunk size available from the workspace.
    //
    INTEGER nb = max((INTEGER)1, min(lwork, lwkopt) / nq);
    //
    INTEGER i = 0;
    INTEGER len = 0;
    INTEGER ldwork = 0;
    if (left) {
        if (notran) {
            for (i = 1; i <= n; i = i + nb) {
                len = min(nb, n - i + 1);
                ldwork = m;
                //
                //              Multiply bottom part of C by Q12.
                //
                Rlacpy("All", n1, len, &c[((n2 + 1) - 1) + (i - 1) * ldc], ldc, work, ldwork);
                Rtrmm("Left", "Lower", "No Transpose", "Non-Unit", n1, len, one, &q[((n2 + 1) - 1) * ldq], ldq, work, ldwork);
                //
                //              Multiply top part of C by Q11.
                //
                Rgemm("No Transpose", "No Transpose", n1, len, n2, one, q, ldq, &c[(i - 1) * ldc], ldc, one, work, ldwork);
                //
                //              Multiply top part of C by Q21.
                //
                Rlacpy("All", n2, len, &c[(i - 1) * ldc], ldc, &work[(n1 + 1) - 1], ldwork);
                Rtrmm("Left", "Upper", "No Transpose", "Non-Unit", n2, len, one, &q[((n1 + 1) - 1)], ldq, &work[(n1 + 1) - 1], ldwork);
                //
                //              Multiply bottom part of C by Q22.
                //
                Rgemm("No Transpose", "No Transpose", n2, len, n1, one, &q[((n1 + 1) - 1) + ((n2 + 1) - 1) * ldq], ldq, &c[((n2 + 1) - 1) + (i - 1) * ldc], ldc, one, &work[(n1 + 1) - 1], ldwork);
                //
                //              Copy everything back.
                //
                Rlacpy("All", m, len, work, ldwork, &c[(i - 1) * ldc], ldc);
            }
        } else {
            for (i = 1; i <= n; i = i + nb) {
                len = min(nb, n - i + 1);
                ldwork = m;
                //
                //              Multiply bottom part of C by Q21**T.
                //
                Rlacpy("All", n2, len, &c[((n1 + 1) - 1) + (i - 1) * ldc], ldc, work, ldwork);
                Rtrmm("Left", "Upper", "Transpose", "Non-Unit", n2, len, one, &q[((n1 + 1) - 1)], ldq, work, ldwork);
                //
                //              Multiply top part of C by Q11**T.
                //
                Rgemm("Transpose", "No Transpose", n2, len, n1, one, q, ldq, &c[(i - 1) * ldc], ldc, one, work, ldwork);
                //
                //              Multiply top part of C by Q12**T.
                //
                Rlacpy("All", n1, len, &c[(i - 1) * ldc], ldc, &work[(n2 + 1) - 1], ldwork);
                Rtrmm("Left", "Lower", "Transpose", "Non-Unit", n1, len, one, &q[((n2 + 1) - 1) * ldq], ldq, &work[(n2 + 1) - 1], ldwork);
                //
                //              Multiply bottom part of C by Q22**T.
                //
                Rgemm("Transpose", "No Transpose", n1, len, n2, one, &q[((n1 + 1) - 1) + ((n2 + 1) - 1) * ldq], ldq, &c[((n1 + 1) - 1) + (i - 1) * ldc], ldc, one, &work[(n2 + 1) - 1], ldwork);
                //
                //              Copy everything back.
                //
                Rlacpy("All", m, len, work, ldwork, &c[(i - 1) * ldc], ldc);
            }
        }
    } else {
        if (notran) {
            for (i = 1; i <= m; i = i + nb) {
                len = min(nb, m - i + 1);
                ldwork = len;
                //
                //              Multiply right part of C by Q21.
                //
                Rlacpy("All", len, n2, &c[(i - 1) + ((n1 + 1) - 1) * ldc], ldc, work, ldwork);
                Rtrmm("Right", "Upper", "No Transpose", "Non-Unit", len, n2, one, &q[((n1 + 1) - 1)], ldq, work, ldwork);
                //
                //              Multiply left part of C by Q11.
                //
                Rgemm("No Transpose", "No Transpose", len, n2, n1, one, &c[(i - 1)], ldc, q, ldq, one, work, ldwork);
                //
                //              Multiply left part of C by Q12.
                //
                Rlacpy("All", len, n1, &c[(i - 1)], ldc, &work[(1 + n2 * ldwork) - 1], ldwork);
                Rtrmm("Right", "Lower", "No Transpose", "Non-Unit", len, n1, one, &q[((n2 + 1) - 1) * ldq], ldq, &work[(1 + n2 * ldwork) - 1], ldwork);
                //
                //              Multiply right part of C by Q22.
                //
                Rgemm("No Transpose", "No Transpose", len, n1, n2, one, &c[(i - 1) + ((n1 + 1) - 1) * ldc], ldc, &q[((n1 + 1) - 1) + ((n2 + 1) - 1) * ldq], ldq, one, &work[(1 + n2 * ldwork) - 1], ldwork);
                //
                //              Copy everything back.
                //
                Rlacpy("All", len, n, work, ldwork, &c[(i - 1)], ldc);
            }
        } else {
            for (i = 1; i <= m; i = i + nb) {
                len = min(nb, m - i + 1);
                ldwork = len;
                //
                //              Multiply right part of C by Q12**T.
                //
                Rlacpy("All", len, n1, &c[(i - 1) + ((n2 + 1) - 1) * ldc], ldc, work, ldwork);
                Rtrmm("Right", "Lower", "Transpose", "Non-Unit", len, n1, one, &q[((n2 + 1) - 1) * ldq], ldq, work, ldwork);
                //
                //              Multiply left part of C by Q11**T.
                //
                Rgemm("No Transpose", "Transpose", len, n1, n2, one, &c[(i - 1)], ldc, q, ldq, one, work, ldwork);
                //
                //              Multiply left part of C by Q21**T.
                //
                Rlacpy("All", len, n2, &c[(i - 1)], ldc, &work[(1 + n1 * ldwork) - 1], ldwork);
                Rtrmm("Right", "Upper", "Transpose", "Non-Unit", len, n2, one, &q[((n1 + 1) - 1)], ldq, &work[(1 + n1 * ldwork) - 1], ldwork);
                //
                //              Multiply right part of C by Q22**T.
                //
                Rgemm("No Transpose", "Transpose", len, n2, n1, one, &c[(i - 1) + ((n2 + 1) - 1) * ldc], ldc, &q[((n1 + 1) - 1) + ((n2 + 1) - 1) * ldq], ldq, one, &work[(1 + n1 * ldwork) - 1], ldwork);
                //
                //              Copy everything back.
                //
                Rlacpy("All", len, n, work, ldwork, &c[(i - 1)], ldc);
            }
        }
    }
    //
    work[1 - 1] = castREAL(lwkopt);
    //
    //     End of Rorm22
    //
}

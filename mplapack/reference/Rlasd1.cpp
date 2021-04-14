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

void Rlasd1(INTEGER const nl, INTEGER const nr, INTEGER const sqre, REAL *d, REAL &alpha, REAL &beta, REAL *u, INTEGER const ldu, REAL *vt, INTEGER const ldvt, INTEGER *idxq, INTEGER *iwork, REAL *work, INTEGER &info) {
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
    //
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    //
    if (nl < 1) {
        info = -1;
    } else if (nr < 1) {
        info = -2;
    } else if ((sqre < 0) || (sqre > 1)) {
        info = -3;
    }
    if (info != 0) {
        Mxerbla("Rlasd1", -info);
        return;
    }
    //
    INTEGER n = nl + nr + 1;
    INTEGER m = n + sqre;
    //
    //     The following values are for bookkeeping purposes only.  They are
    //     integer pointers which indicate the portion of the workspace
    //     used by a particular array in Rlasd2 and Rlasd3.
    //
    INTEGER ldu2 = n;
    INTEGER ldvt2 = m;
    //
    INTEGER iz = 1;
    INTEGER isigma = iz + m;
    INTEGER iu2 = isigma + n;
    INTEGER ivt2 = iu2 + ldu2 * n;
    INTEGER iq = ivt2 + ldvt2 * m;
    //
    INTEGER idx = 1;
    INTEGER idxc = idx + n;
    INTEGER coltyp = idxc + n;
    INTEGER idxp = coltyp + n;
    //
    //     Scale.
    //
    REAL orgnrm = max(abs(alpha), abs(beta));
    const REAL zero = 0.0;
    d[(nl + 1) - 1] = zero;
    INTEGER i = 0;
    for (i = 1; i <= n; i = i + 1) {
        if (abs(d[i - 1]) > orgnrm) {
            orgnrm = abs(d[i - 1]);
        }
    }
    const REAL one = 1.0;
    Rlascl("G", 0, 0, orgnrm, one, n, 1, d, n, info);
    alpha = alpha / orgnrm;
    beta = beta / orgnrm;
    //
    //     Deflate singular values.
    //
    INTEGER k = 0;
    Rlasd2(nl, nr, sqre, k, d, &work[iz - 1], alpha, beta, u, ldu, vt, ldvt, &work[isigma - 1], &work[iu2 - 1], ldu2, &work[ivt2 - 1], ldvt2, iwork[idxp - 1], iwork[idx - 1], iwork[idxc - 1], idxq, iwork[coltyp - 1], info);
    //
    //     Solve Secular Equation and update singular vectors.
    //
    INTEGER ldq = k;
    Rlasd3(nl, nr, sqre, k, d, &work[iq - 1], ldq, &work[isigma - 1], u, ldu, &work[iu2 - 1], ldu2, vt, ldvt, &work[ivt2 - 1], ldvt2, iwork[idxc - 1], iwork[coltyp - 1], &work[iz - 1], info);
    //
    //     Report the convergence failure.
    //
    if (info != 0) {
        return;
    }
    //
    //     Unscale.
    //
    Rlascl("G", 0, 0, one, orgnrm, n, 1, d, n, info);
    //
    //     Prepare the IDXQ sorting permutation.
    //
    INTEGER n1 = k;
    INTEGER n2 = n - k;
    Rlamrg(n1, n2, d, 1, -1, idxq);
    //
    //     End of Rlasd1
    //
}

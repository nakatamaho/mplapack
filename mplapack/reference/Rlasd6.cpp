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

void Rlasd6(INTEGER const &icompq, INTEGER const &nl, INTEGER const &nr, INTEGER const &sqre, REAL *d, REAL *vf, REAL *vl, REAL &alpha, REAL &beta, INTEGER *idxq, INTEGER *perm, INTEGER const &givptr, arr_cref<INTEGER, 2> givcol, INTEGER const &ldgcol, REAL *givnum, INTEGER const &ldgnum, REAL *poles, REAL *difl, REAL *difr, REAL *z, INTEGER const &k, REAL const &c, REAL const &s, REAL *work, INTEGER *iwork, INTEGER &info) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    INTEGER n = nl + nr + 1;
    INTEGER m = n + sqre;
    //
    if ((icompq < 0) || (icompq > 1)) {
        info = -1;
    } else if (nl < 1) {
        info = -2;
    } else if (nr < 1) {
        info = -3;
    } else if ((sqre < 0) || (sqre > 1)) {
        info = -4;
    } else if (ldgcol < n) {
        info = -14;
    } else if (ldgnum < n) {
        info = -16;
    }
    if (info != 0) {
        Mxerbla("Rlasd6", -info);
        return;
    }
    //
    //     The following values are for bookkeeping purposes only.  They are
    //     INTEGEReger poINTEGERers which indicate the portion of the workspace
    //     used by a particular array in Rlasd7 and Rlasd8.
    //
    INTEGER isigma = 1;
    INTEGER iw = isigma + n;
    INTEGER ivfw = iw + m;
    INTEGER ivlw = ivfw + m;
    //
    INTEGER idx = 1;
    INTEGER idxc = idx + n;
    INTEGER idxp = idxc + n;
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
    //     Sort and Deflate singular values.
    //
    Rlasd7(icompq, nl, nr, sqre, k, d, z, work[iw - 1], vf, work[ivfw - 1], vl, work[ivlw - 1], alpha, beta, work[isigma - 1], iwork[idx - 1], iwork[idxp - 1], idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info);
    //
    //     Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.
    //
    Rlasd8(icompq, k, d, z, vf, vl, difl, difr, ldgnum, work[isigma - 1], work[iw - 1], info);
    //
    //     Report the possible convergence failure.
    //
    if (info != 0) {
        return;
    }
    //
    //     Save the poles if ICOMPQ = 1.
    //
    if (icompq == 1) {
        Rcopy(k, d, 1, poles[(1 - 1)], 1);
        Rcopy(k, work[isigma - 1], 1, poles[(2 - 1) * ldpoles], 1);
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
    //     End of Rlasd6
    //
}

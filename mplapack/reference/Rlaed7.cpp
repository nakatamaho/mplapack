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

void Rlaed7(INTEGER const icompq, INTEGER const n, INTEGER const qsiz, INTEGER const tlvls, INTEGER const curlvl, INTEGER const curpbm, REAL *d, REAL *q, INTEGER const ldq, INTEGER *indxq, REAL rho, INTEGER const cutpnt, REAL *qstore, INTEGER *qptr, INTEGER *prmptr, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, REAL *givnum, REAL *work, INTEGER *iwork, INTEGER &info) {
    INTEGER ldq2 = 0;
    INTEGER iz = 0;
    INTEGER idlmda = 0;
    INTEGER iw = 0;
    INTEGER iq2 = 0;
    INTEGER is = 0;
    INTEGER indx = 0;
    INTEGER indxc = 0;
    INTEGER coltyp = 0;
    INTEGER indxp = 0;
    INTEGER ptr = 0;
    INTEGER i = 0;
    INTEGER curr = 0;
    INTEGER k = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    //
    if (icompq < 0 || icompq > 1) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (icompq == 1 && qsiz < n) {
        info = -3;
    } else if (ldq < max((INTEGER)1, n)) {
        info = -9;
    } else if (min(1, n) > cutpnt || n < cutpnt) {
        info = -12;
    }
    if (info != 0) {
        Mxerbla("Rlaed7", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     The following values are for bookkeeping purposes only.  They are
    //     integer pointers which indicate the portion of the workspace
    //     used by a particular array in DLAED8 and DLAED9.
    //
    if (icompq == 1) {
        ldq2 = qsiz;
    } else {
        ldq2 = n;
    }
    //
    iz = 1;
    idlmda = iz + n;
    iw = idlmda + n;
    iq2 = iw + n;
    is = iq2 + n * ldq2;
    //
    indx = 1;
    indxc = indx + n;
    coltyp = indxc + n;
    indxp = coltyp + n;
    //
    //     Form the z-vector which consists of the last row of Q_1 and the
    //     first row of Q_2.
    //
    ptr = 1 + pow(2, tlvls);
    for (i = 1; i <= curlvl - 1; i = i + 1) {
        ptr += pow(2, (tlvls - i));
    }
    curr = ptr + curpbm;
    Rlaeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, qstore, qptr, &work[iz - 1], &work[(iz + n) - 1], info);
    //
    //     When solving the final problem, we no longer need the stored data,
    //     so we will overwrite the data from this level onto the previously
    //     used storage space.
    //
    if (curlvl == tlvls) {
        qptr[curr - 1] = 1;
        prmptr[curr - 1] = 1;
        givptr[curr - 1] = 1;
    }
    //
    //     Sort and Deflate eigenvalues.
    //
    Rlaed8(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, &work[iz - 1], &work[idlmda - 1], &work[iq2 - 1], ldq2, &work[iw - 1], &perm[prmptr[curr - 1] - 1], givptr[(curr + 1) - 1], &givcol[(givptr[curr - 1] - 1) * 2], &givnum[(givptr[curr - 1] - 1) * 2], &iwork[indxp - 1], &iwork[indx - 1], info);
    prmptr[(curr + 1) - 1] = prmptr[curr - 1] + n;
    givptr[(curr + 1) - 1] += givptr[curr - 1];
    //
    //     Solve Secular Equation.
    //
    if (k != 0) {
        Rlaed9(k, 1, k, n, d, &work[is - 1], k, rho, &work[idlmda - 1], &work[iw - 1], &qstore[qptr[curr - 1] - 1], k, info);
        if (info != 0) {
            goto statement_30;
        }
        if (icompq == 1) {
            Rgemm("N", "N", qsiz, k, k, one, &work[iq2 - 1], ldq2, &qstore[qptr[curr - 1] - 1], k, zero, q, ldq);
        }
        qptr[(curr + 1) - 1] = qptr[curr - 1] + pow2(k);
        //
        //     Prepare the INDXQ sorting permutation.
        //
        n1 = k;
        n2 = n - k;
        Rlamrg(n1, n2, d, 1, -1, indxq);
    } else {
        qptr[(curr + 1) - 1] = qptr[curr - 1];
        for (i = 1; i <= n; i = i + 1) {
            indxq[i - 1] = i;
        }
    }
//
statement_30:;
    //
    //     End of Rlaed7
    //
}

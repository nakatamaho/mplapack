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

void Claed7(INTEGER const &n, INTEGER const &cutpnt, INTEGER const &qsiz, INTEGER const &tlvls, INTEGER const &curlvl, INTEGER const &curpbm, REAL *d, COMPLEX *q, INTEGER const &ldq, REAL const &rho, arr_ref<INTEGER> indxq, REAL *qstore, arr_ref<INTEGER> qptr, arr_ref<INTEGER> prmptr, INTEGER *perm, arr_ref<INTEGER> givptr, arr_cref<INTEGER, 2> givcol, REAL *givnum, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER &info) {
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
    //     IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
    //        INFO = -1
    //     ELSE IF( N.LT.0 ) THEN
    if (n < 0) {
        info = -1;
    } else if (min(1, n) > cutpnt || n < cutpnt) {
        info = -2;
    } else if (qsiz < n) {
        info = -3;
    } else if (ldq < max((INTEGER)1, n)) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Claed7", -info);
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
    //     INTEGEReger poINTEGERers which indicate the portion of the workspace
    //     used by a particular array in Rlaed2 and SLAED3.
    //
    INTEGER iz = 1;
    INTEGER idlmda = iz + n;
    INTEGER iw = idlmda + n;
    INTEGER iq = iw + n;
    //
    INTEGER indx = 1;
    INTEGER indxc = indx + n;
    INTEGER coltyp = indxc + n;
    INTEGER indxp = coltyp + n;
    //
    //     Form the z-vector which consists of the last row of Q_1 and the
    //     first row of Q_2.
    //
    INTEGER ptr = 1 + pow(2, tlvls);
    INTEGER i = 0;
    for (i = 1; i <= curlvl - 1; i = i + 1) {
        ptr += pow(2, (tlvls - i));
    }
    INTEGER curr = ptr + curpbm;
    Rlaeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, qstore, qptr, rwork[iz - 1], rwork[(iz + n) - 1], info);
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
    INTEGER k = 0;
    Claed8(k, n, qsiz, q, ldq, d, rho, cutpnt, rwork[iz - 1], rwork[idlmda - 1], work, qsiz, rwork[iw - 1], iwork[indxp - 1], iwork[indx - 1], indxq, perm[prmptr[curr - 1] - 1], givptr[(curr + 1) - 1], givcol[(givptr[curr - 1] - 1) * ldgivcol], givnum[(givptr[curr - 1] - 1) * ldgivnum], info);
    prmptr[(curr + 1) - 1] = prmptr[curr - 1] + n;
    givptr[(curr + 1) - 1] += givptr[curr - 1];
    //
    //     Solve Secular Equation.
    //
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    if (k != 0) {
        Rlaed9(k, 1, k, n, d, rwork[iq - 1], k, rho, rwork[idlmda - 1], rwork[iw - 1], qstore[qptr[curr - 1] - 1], k, info);
        Clacrm(qsiz, k, work, qsiz, qstore[qptr[curr - 1] - 1], k, q, ldq, rwork[iq - 1]);
        qptr[(curr + 1) - 1] = qptr[curr - 1] + pow2(k);
        if (info != 0) {
            return;
        }
        //
        //     Prepare the INDXQ sorting premutation.
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
    //     End of Claed7
    //
}

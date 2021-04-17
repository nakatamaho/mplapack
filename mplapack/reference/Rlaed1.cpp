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

void Rlaed1(INTEGER const n, REAL *d, REAL *q, INTEGER const ldq, INTEGER *indxq, REAL rho, INTEGER const cutpnt, REAL *work, INTEGER *iwork, INTEGER &info) {
    INTEGER iz = 0;
    INTEGER idlmda = 0;
    INTEGER iw = 0;
    INTEGER iq2 = 0;
    INTEGER indx = 0;
    INTEGER indxc = 0;
    INTEGER coltyp = 0;
    INTEGER indxp = 0;
    INTEGER zpp1 = 0;
    INTEGER k = 0;
    INTEGER is = 0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER i = 0;
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
    if (n < 0) {
        info = -1;
    } else if (ldq < max((INTEGER)1, n)) {
        info = -4;
    } else if (min((INTEGER)1, n / 2) > cutpnt || (n / 2) < cutpnt) {
        info = -7;
    }
    if (info != 0) {
        Mxerbla("Rlaed1", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     The following values are integer pointers which indicate
    //     the portion of the workspace
    //     used by a particular array in Rlaed2 and Rlaed3.
    //
    iz = 1;
    idlmda = iz + n;
    iw = idlmda + n;
    iq2 = iw + n;
    //
    indx = 1;
    indxc = indx + n;
    coltyp = indxc + n;
    indxp = coltyp + n;
    //
    //     Form the z-vector which consists of the last row of Q_1 and the
    //     first row of Q_2.
    //
    Rcopy(cutpnt, &q[(cutpnt - 1)], ldq, &work[iz - 1], 1);
    zpp1 = cutpnt + 1;
    Rcopy(n - cutpnt, &q[(zpp1 - 1) + (zpp1 - 1) * ldq], ldq, &work[(iz + cutpnt) - 1], 1);
    //
    //     Deflate eigenvalues.
    //
    Rlaed2(k, n, cutpnt, d, q, ldq, indxq, rho, &work[iz - 1], &work[idlmda - 1], &work[iw - 1], &work[iq2 - 1], &iwork[indx - 1], &iwork[indxc - 1], &iwork[indxp - 1], &iwork[coltyp - 1], info);
    //
    if (info != 0) {
        goto statement_20;
    }
    //
    //     Solve Secular Equation.
    //
    if (k != 0) {
        is = (iwork[coltyp - 1] + iwork[(coltyp + 1) - 1]) * cutpnt + (iwork[(coltyp + 1) - 1] + iwork[(coltyp + 2) - 1]) * (n - cutpnt) + iq2;
        Rlaed3(k, n, cutpnt, d, q, ldq, rho, &work[idlmda - 1], &work[iq2 - 1], &iwork[indxc - 1], &iwork[coltyp - 1], &work[iw - 1], &work[is - 1], info);
        if (info != 0) {
            goto statement_20;
        }
        //
        //     Prepare the INDXQ sorting permutation.
        //
        n1 = k;
        n2 = n - k;
        Rlamrg(n1, n2, d, 1, -1, indxq);
    } else {
        for (i = 1; i <= n; i = i + 1) {
            indxq[i - 1] = i;
        }
    }
//
statement_20:;
    //
    //     End of Rlaed1
    //
}

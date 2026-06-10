/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DLAEDA.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Rlaeda(INTEGER const n, INTEGER const tlvls, INTEGER const curlvl, INTEGER const curpbm, INTEGER *prmptr, INTEGER *perm, INTEGER *givptr, INTEGER *givcol, REAL *givnum, REAL *q, INTEGER *qptr, REAL *z, REAL *ztemp, INTEGER &info) {
    INTEGER ldgivcol = 2;
    INTEGER ldgivnum = 2;
    //
    // Test the input parameters.
    //
    info = 0;
    //
    if (n < 0) {
        info = -1;
    }
    if (info != 0) {
        Mxerbla("Rlaeda", -info);
        return;
    }
    //
    // Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    // Determine location of first number in second half.
    //
    INTEGER mid = n / 2 + 1;
    //
    // Gather last/first rows of appropriate eigenblocks into center of Z
    //
    INTEGER ptr = 1;
    //
    // Determine location of lowest level subproblem in the full storage
    // scheme
    //
    INTEGER curr = ptr + curpbm * (INTEGER(1) << (curlvl)) + (INTEGER(1) << ((curlvl - 1))) - 1;
    //
    // Determine size of these matrices.  We add HALF to the value of
    // the SQRT in case the machine underestimates one of these square
    // roots.
    //
    const REAL half = 0.5;
    INTEGER bsiz1 = castINTEGER(half + sqrt(castREAL(qptr[(curr + 1) - 1] - qptr[curr - 1])));
    INTEGER bsiz2 = castINTEGER(half + sqrt(castREAL(qptr[(curr + 2) - 1] - qptr[(curr + 1) - 1])));
    INTEGER k = 0;
    const REAL zero = 0.0;
    for (k = 1; k <= mid - bsiz1 - 1; k = k + 1) {
        z[k - 1] = zero;
    }
    Rcopy(bsiz1, &q[(qptr[curr - 1] + bsiz1 - 1) - 1], bsiz1, &z[(mid - bsiz1) - 1], 1);
    Rcopy(bsiz2, &q[qptr[(curr + 1) - 1] - 1], bsiz2, &z[mid - 1], 1);
    for (k = mid + bsiz2; k <= n; k = k + 1) {
        z[k - 1] = zero;
    }
    //
    // Loop through remaining levels 1 -> CURLVL applying the Givens
    // rotations and permutation and then multiplying the center matrices
    // against the current Z.
    //
    ptr = (INTEGER(1) << (tlvls)) + 1;
    INTEGER psiz1 = 0;
    INTEGER psiz2 = 0;
    INTEGER zptr1 = 0;
    INTEGER i = 0;
    const REAL one = 1.0;
    for (k = 1; k <= curlvl - 1; k = k + 1) {
        curr = ptr + curpbm * (INTEGER(1) << ((curlvl - k))) + (INTEGER(1) << ((curlvl - k - 1))) - 1;
        psiz1 = prmptr[(curr + 1) - 1] - prmptr[curr - 1];
        psiz2 = prmptr[(curr + 2) - 1] - prmptr[(curr + 1) - 1];
        zptr1 = mid - psiz1;
        //
        // Apply Givens at CURR and CURR+1
        //
        for (i = givptr[curr - 1]; i <= givptr[(curr + 1) - 1] - 1; i = i + 1) {
            Rrot(1, &z[(zptr1 + givcol[(i - 1) * ldgivcol] - 1) - 1], 1, &z[(zptr1 + givcol[(2 - 1) + (i - 1) * ldgivcol] - 1) - 1], 1, givnum[(i - 1) * ldgivnum], givnum[(2 - 1) + (i - 1) * ldgivnum]);
        }
        for (i = givptr[(curr + 1) - 1]; i <= givptr[(curr + 2) - 1] - 1; i = i + 1) {
            Rrot(1, &z[(mid - 1 + givcol[(i - 1) * ldgivcol]) - 1], 1, &z[(mid - 1 + givcol[(2 - 1) + (i - 1) * ldgivcol]) - 1], 1, givnum[(i - 1) * ldgivnum], givnum[(2 - 1) + (i - 1) * ldgivnum]);
        }
        psiz1 = prmptr[(curr + 1) - 1] - prmptr[curr - 1];
        psiz2 = prmptr[(curr + 2) - 1] - prmptr[(curr + 1) - 1];
        for (i = 0; i <= psiz1 - 1; i = i + 1) {
            ztemp[(i + 1) - 1] = z[(zptr1 + perm[(prmptr[curr - 1] + i) - 1] - 1) - 1];
        }
        for (i = 0; i <= psiz2 - 1; i = i + 1) {
            ztemp[(psiz1 + i + 1) - 1] = z[(mid + perm[(prmptr[(curr + 1) - 1] + i) - 1] - 1) - 1];
        }
        //
        // Multiply Blocks at CURR and CURR+1
        //
        // Determine size of these matrices.  We add HALF to the value of
        // the SQRT in case the machine underestimates one of these
        // square roots.
        //
        bsiz1 = castINTEGER(half + sqrt(castREAL(qptr[(curr + 1) - 1] - qptr[curr - 1])));
        bsiz2 = castINTEGER(half + sqrt(castREAL(qptr[(curr + 2) - 1] - qptr[(curr + 1) - 1])));
        if (bsiz1 > 0) {
            Rgemv("T", bsiz1, bsiz1, one, &q[qptr[curr - 1] - 1], bsiz1, &ztemp[1 - 1], 1, zero, &z[zptr1 - 1], 1);
        }
        Rcopy(psiz1 - bsiz1, &ztemp[(bsiz1 + 1) - 1], 1, &z[(zptr1 + bsiz1) - 1], 1);
        if (bsiz2 > 0) {
            Rgemv("T", bsiz2, bsiz2, one, &q[qptr[(curr + 1) - 1] - 1], bsiz2, &ztemp[(psiz1 + 1) - 1], 1, zero, &z[mid - 1], 1);
        }
        Rcopy(psiz2 - bsiz2, &ztemp[(psiz1 + bsiz2 + 1) - 1], 1, &z[(mid + bsiz2) - 1], 1);
        //
        ptr += (INTEGER(1) << ((tlvls - k)));
    }
    //
    // End of Rlaeda
    //
}

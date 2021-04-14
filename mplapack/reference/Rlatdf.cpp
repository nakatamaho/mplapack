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

void Rlatdf(INTEGER const ijob, INTEGER const n, REAL *z, INTEGER const ldz, REAL *rhs, REAL const rdsum, REAL const rRscal, INTEGER *ipiv, INTEGER *jpiv) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    const REAL one = 1.0;
    REAL pmone = 0.0;
    INTEGER j = 0;
    REAL bp = 0.0;
    REAL bm = 0.0;
    REAL splus = 0.0;
    REAL sminu = 0.0;
    REAL temp = 0.0;
    const INTEGER maxdim = 8;
    arr_1d<maxdim, REAL> xp(fill0);
    const REAL zero = 0.0;
    INTEGER i = 0;
    INTEGER k = 0;
    arr_1d<4 * maxdim, REAL> work(fill0);
    arr_1d<maxdim, int> iwork(fill0);
    INTEGER info = 0;
    arr_1d<maxdim, REAL> xm(fill0);
    if (ijob != 2) {
        //
        //        Apply permutations IPIV to RHS
        //
        Rlaswp(1, rhs, ldz, 1, n - 1, ipiv, 1);
        //
        //        Solve for L-part choosing RHS either to +1 or -1.
        //
        pmone = -one;
        //
        for (j = 1; j <= n - 1; j = j + 1) {
            bp = rhs[j - 1] + one;
            bm = rhs[j - 1] - one;
            splus = one;
            //
            //           Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and
            //           SMIN computed more efficiently than in BSOLVE [1].
            //
            splus += Rdot(n - j, &z[((j + 1) - 1) + (j - 1) * ldz], 1, &z[((j + 1) - 1) + (j - 1) * ldz], 1);
            sminu = Rdot(n - j, &z[((j + 1) - 1) + (j - 1) * ldz], 1, rhs[(j + 1) - 1], 1);
            splus = splus * rhs[j - 1];
            if (splus > sminu) {
                rhs[j - 1] = bp;
            } else if (sminu > splus) {
                rhs[j - 1] = bm;
            } else {
                //
                //              In this case the updating sums are equal and we can
                //              choose RHS(J) +1 or -1. The first time this happens
                //              we choose -1, thereafter +1. This is a simple way to
                //              get good estimates of matrices like Byers well-known
                //              example (see [1]). (Not done in BSOLVE.)
                //
                rhs[j - 1] += pmone;
                pmone = one;
            }
            //
            //           Compute the remaining r.h.s.
            //
            temp = -rhs[j - 1];
            Raxpy(n - j, temp, &z[((j + 1) - 1) + (j - 1) * ldz], 1, rhs[(j + 1) - 1], 1);
            //
        }
        //
        //        Solve for U-part, look-ahead for RHS(N) = +-1. This is not done
        //        in BSOLVE and will hopefully give us a better estimate because
        //        any ill-conditioning of the original matrix is transferred to U
        //        and not to L. U(N, N) is an approximation to sigma_min(LU).
        //
        Rcopy(n - 1, rhs, 1, xp, 1);
        xp[n - 1] = rhs[n - 1] + one;
        rhs[n - 1] = rhs[n - 1] - one;
        splus = zero;
        sminu = zero;
        for (i = n; i >= 1; i = i - 1) {
            temp = one / z[(i - 1) + (i - 1) * ldz];
            xp[i - 1] = xp[i - 1] * temp;
            rhs[i - 1] = rhs[i - 1] * temp;
            for (k = i + 1; k <= n; k = k + 1) {
                xp[i - 1] = xp[i - 1] - xp[k - 1] * (z[(i - 1) + (k - 1) * ldz] * temp);
                rhs[i - 1] = rhs[i - 1] - rhs[k - 1] * (z[(i - 1) + (k - 1) * ldz] * temp);
            }
            splus += abs(xp[i - 1]);
            sminu += abs(rhs[i - 1]);
        }
        if (splus > sminu) {
            Rcopy(n, xp, 1, rhs, 1);
        }
        //
        //        Apply the permutations JPIV to the computed solution (RHS)
        //
        Rlaswp(1, rhs, ldz, 1, n - 1, jpiv, -1);
        //
        //        Compute the sum of squares
        //
        Rlassq(n, rhs, 1, rRscal, rdsum);
        //
    } else {
        //
        //        IJOB = 2, Compute approximate nullvector XM of Z
        //
        Rgecon("I", n, z, ldz, one, temp, work, iwork, info);
        Rcopy(n, &work[(n + 1) - 1], 1, xm, 1);
        //
        //        Compute RHS
        //
        Rlaswp(1, xm, ldz, 1, n - 1, ipiv, -1);
        temp = one / sqrt(Rdot(n, xm, 1, xm, 1));
        Rscal(n, temp, xm, 1);
        Rcopy(n, xm, 1, xp, 1);
        Raxpy(n, one, rhs, 1, xp, 1);
        Raxpy(n, -one, xm, 1, rhs, 1);
        Rgesc2(n, z, ldz, rhs, ipiv, jpiv, temp);
        Rgesc2(n, z, ldz, xp, ipiv, jpiv, temp);
        if (Rasum(n, xp, 1) > Rasum(n, rhs, 1)) {
            Rcopy(n, xp, 1, rhs, 1);
        }
        //
        //        Compute the sum of squares
        //
        Rlassq(n, rhs, 1, rRscal, rdsum);
        //
    }
    //
    //     End of Rlatdf
    //
}

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

void Clatdf(INTEGER const ijob, INTEGER const n, COMPLEX *z, INTEGER const ldz, COMPLEX *rhs, REAL const rdsum, REAL const rRscal, INTEGER *ipiv, INTEGER *jpiv) {
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
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    COMPLEX pmone = 0.0;
    INTEGER j = 0;
    COMPLEX bp = 0.0;
    COMPLEX bm = 0.0;
    const REAL one = 1.0;
    REAL splus = 0.0;
    REAL sminu = 0.0;
    COMPLEX temp = 0.0;
    const INTEGER maxdim = 2;
    arr_1d<4 * maxdim, COMPLEX> work(fill0);
    const REAL zero = 0.0;
    INTEGER i = 0;
    INTEGER k = 0;
    if (ijob != 2) {
        //
        //        Apply permutations IPIV to RHS
        //
        Claswp(1, rhs, ldz, 1, n - 1, ipiv, 1);
        //
        //        Solve for L-part choosing RHS either to +1 or -1.
        //
        pmone = -cone;
        for (j = 1; j <= n - 1; j = j + 1) {
            bp = rhs[j - 1] + cone;
            bm = rhs[j - 1] - cone;
            splus = one;
            //
            //           Lockahead for L- part RHS(1:N-1) = +-1
            //           SPLUS and SMIN computed more efficiently than in BSOLVE[1].
            //
            splus += Cdotc(n - j, &z[((j + 1) - 1) + (j - 1) * ldz], 1, &z[((j + 1) - 1) + (j - 1) * ldz], 1).real();
            sminu = Cdotc(n - j, &z[((j + 1) - 1) + (j - 1) * ldz], 1, rhs[(j + 1) - 1], 1).real();
            splus = splus * rhs[j - 1].real();
            if (splus > sminu) {
                rhs[j - 1] = bp;
            } else if (sminu > splus) {
                rhs[j - 1] = bm;
            } else {
                //
                //              In this case the updating sums are equal and we can
                //              choose RHS(J) +1 or -1. The first time this happens we
                //              choose -1, thereafter +1. This is a simple way to get
                //              good estimates of matrices like Byers well-known example
                //              (see [1]). (Not done in BSOLVE.)
                //
                rhs[j - 1] += pmone;
                pmone = cone;
            }
            //
            //           Compute the remaining r.h.s.
            //
            temp = -rhs[j - 1];
            Caxpy(n - j, temp, &z[((j + 1) - 1) + (j - 1) * ldz], 1, rhs[(j + 1) - 1], 1);
        }
        //
        //        Solve for U- part, lockahead for RHS(N) = +-1. This is not done
        //        In BSOLVE and will hopefully give us a better estimate because
        //        any ill-conditioning of the original matrix is transferred to U
        //        and not to L. U(N, N) is an approximation to sigma_min(LU).
        //
        Ccopy(n - 1, rhs, 1, work, 1);
        work[n - 1] = rhs[n - 1] + cone;
        rhs[n - 1] = rhs[n - 1] - cone;
        splus = zero;
        sminu = zero;
        for (i = n; i >= 1; i = i - 1) {
            temp = cone / z[(i - 1) + (i - 1) * ldz];
            work[i - 1] = work[i - 1] * temp;
            rhs[i - 1] = rhs[i - 1] * temp;
            for (k = i + 1; k <= n; k = k + 1) {
                work[i - 1] = work[i - 1] - work[k - 1] * (z[(i - 1) + (k - 1) * ldz] * temp);
                rhs[i - 1] = rhs[i - 1] - rhs[k - 1] * (z[(i - 1) + (k - 1) * ldz] * temp);
            }
            splus += abs(work[i - 1]);
            sminu += abs(rhs[i - 1]);
        }
        if (splus > sminu) {
            Ccopy(n, work, 1, rhs, 1);
        }
        //
        //        Apply the permutations JPIV to the computed solution (RHS)
        //
        Claswp(1, rhs, ldz, 1, n - 1, jpiv, -1);
        //
        //        Compute the sum of squares
        //
        Classq(n, rhs, 1, rRscal, rdsum);
        return;
    }
    //
    //     ENTRY IJOB = 2
    //
    //     Compute approximate nullvector XM of Z
    //
    REAL rtemp = 0.0;
    arr_1d<maxdim, REAL> rwork(fill0);
    INTEGER info = 0;
    Cgecon("I", n, z, ldz, one, rtemp, work, rwork, info);
    arr_1d<maxdim, COMPLEX> xm(fill0);
    Ccopy(n, &work[(n + 1) - 1], 1, xm, 1);
    //
    //     Compute RHS
    //
    Claswp(1, xm, ldz, 1, n - 1, ipiv, -1);
    temp = cone / sqrt(Cdotc(n, xm, 1, xm, 1));
    Cscal(n, temp, xm, 1);
    arr_1d<maxdim, COMPLEX> xp(fill0);
    Ccopy(n, xm, 1, xp, 1);
    Caxpy(n, cone, rhs, 1, xp, 1);
    Caxpy(n, -cone, xm, 1, rhs, 1);
    REAL scale = 0.0;
    Cgesc2(n, z, ldz, rhs, ipiv, jpiv, scale);
    Cgesc2(n, z, ldz, xp, ipiv, jpiv, scale);
    if (RCasum(n, xp, 1) > RCasum(n, rhs, 1)) {
        Ccopy(n, xp, 1, rhs, 1);
    }
    //
    //     Compute the sum of squares
    //
    Classq(n, rhs, 1, rRscal, rdsum);
    //
    //     End of Clatdf
    //
}

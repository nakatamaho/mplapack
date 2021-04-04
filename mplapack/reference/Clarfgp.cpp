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

void Clarfgp(INTEGER const &n, COMPLEX &alpha, COMPLEX *x, INTEGER const &incx, COMPLEX &tau) {
    const REAL zero = 0.0;
    REAL xnorm = 0.0;
    REAL alphr = 0.0;
    REAL alphi = 0.0;
    const REAL two = 2.0e+0;
    INTEGER j = 0;
    const REAL one = 1.0;
    REAL beta = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    INTEGER knt = 0;
    COMPLEX savealpha = 0.0;
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    if (n <= 0) {
        tau = zero;
        return;
    }
    //
    xnorm = RCnrm2[((n - 1) - 1) + (x - 1) * ldRCnrm2];
    alphr = alpha.real();
    alphi = alpha.imag();
    //
    if (xnorm == zero) {
        //
        //        H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0.
        //
        if (alphi == zero) {
            if (alphr >= zero) {
                //              When TAU.eq.ZERO, the vector is special-cased to be
                //              all zeros in the application routines.  We do not need
                //              to clear it.
                tau = zero;
            } else {
                //              However, the application routines rely on explicit
                //              zero checks when TAU.ne.ZERO, and we must clear X.
                tau = two;
                for (j = 1; j <= n - 1; j = j + 1) {
                    x[(1 + (j - 1) * incx) - 1] = zero;
                }
                alpha = -alpha;
            }
        } else {
            //           Only "reflecting" the diagonal entry to be real and non-negative.
            xnorm = Rlapy2[(alphr - 1) + (alphi - 1) * ldRlapy2];
            tau = COMPLEX(one - alphr / xnorm, -alphi / xnorm);
            for (j = 1; j <= n - 1; j = j + 1) {
                x[(1 + (j - 1) * incx) - 1] = zero;
            }
            alpha = xnorm;
        }
    } else {
        //
        //        general case
        //
        beta = sign[(Rlapy3[(alphr - 1) + (alphi - 1) * ldRlapy3] - 1) + (alphr - 1) * ldsign];
        smlnum = dlamch("S") / dlamch("E");
        bignum = one / smlnum;
        //
        knt = 0;
        if (abs(beta) < smlnum) {
        //
        //           XNORM, BETA may be inaccurate; scale X and recompute them
        //
        statement_10:
            knt++;
            CRscal(n - 1, bignum, x, incx);
            beta = beta * bignum;
            alphi = alphi * bignum;
            alphr = alphr * bignum;
            if ((abs(beta) < smlnum) && (knt < 20)) {
                goto statement_10;
            }
            //
            //           New BETA is at most 1, at least SMLNUM
            //
            xnorm = RCnrm2[((n - 1) - 1) + (x - 1) * ldRCnrm2];
            alpha = COMPLEX(alphr, alphi);
            beta = sign[(Rlapy3[(alphr - 1) + (alphi - 1) * ldRlapy3] - 1) + (alphr - 1) * ldsign];
        }
        savealpha = alpha;
        alpha += beta;
        if (beta < zero) {
            beta = -beta;
            tau = -alpha / beta;
        } else {
            alphr = alphi * (alphi / alpha.real());
            alphr += xnorm * (xnorm / alpha.real());
            tau = COMPLEX(alphr / beta, -alphi / beta);
            alpha = COMPLEX(-alphr, alphi);
        }
        alpha = Cladiv[(COMPLEX(one) - 1) + (alpha - 1) * ldCladiv];
        //
        if (abs(tau) <= smlnum) {
            //
            //           In the case where the computed TAU ends up being a denormalized number,
            //           it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            //           to ZERO (or TWO or whatever makes a nonnegative real number for BETA).
            //
            //           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            //           (Thanks Pat. Thanks MathWorks.)
            //
            alphr = savealpha.real();
            alphi = savealpha.imag();
            if (alphi == zero) {
                if (alphr >= zero) {
                    tau = zero;
                } else {
                    tau = two;
                    for (j = 1; j <= n - 1; j = j + 1) {
                        x[(1 + (j - 1) * incx) - 1] = zero;
                    }
                    beta = -savealpha;
                }
            } else {
                xnorm = Rlapy2[(alphr - 1) + (alphi - 1) * ldRlapy2];
                tau = COMPLEX(one - alphr / xnorm, -alphi / xnorm);
                for (j = 1; j <= n - 1; j = j + 1) {
                    x[(1 + (j - 1) * incx) - 1] = zero;
                }
                beta = xnorm;
            }
            //
        } else {
            //
            //           This is the general case.
            //
            Cscal(n - 1, alpha, x, incx);
            //
        }
        //
        //        If BETA is subnormal, it may lose relative accuracy
        //
        for (j = 1; j <= knt; j = j + 1) {
            beta = beta * smlnum;
        }
        alpha = beta;
    }
    //
    //     End of Clarfgp
    //
}

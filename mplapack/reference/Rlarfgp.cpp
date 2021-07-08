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

void Rlarfgp(INTEGER const n, REAL &alpha, REAL *x, INTEGER const incx, REAL &tau) {
    const REAL zero = 0.0;
    REAL xnorm = 0.0;
    const REAL two = 2.0e+0;
    INTEGER j = 0;
    REAL beta = 0.0;
    REAL smlnum = 0.0;
    INTEGER knt = 0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL savealpha = 0.0;
    //
    if (n <= 0) {
        tau = zero;
        return;
    }
    //
    xnorm = Rnrm2(n - 1, x, incx);
    //
    if (xnorm == zero) {
        //
        //        H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0
        //
        if (alpha >= zero) {
            //           When TAU.eq.ZERO, the vector is special-cased to be
            //           all zeros in the application routines.  We do not need
            //           to clear it.
            tau = zero;
        } else {
            //           However, the application routines rely on explicit
            //           zero checks when TAU.ne.ZERO, and we must clear X.
            tau = two;
            for (j = 1; j <= n - 1; j = j + 1) {
                x[(1 + (j - 1) * incx) - 1] = 0.0;
            }
            alpha = -alpha;
        }
    } else {
        //
        //        general case
        //
        beta = sign(Rlapy2(alpha, xnorm), alpha);
        smlnum = Rlamch("S") / Rlamch("E");
        knt = 0;
        if (abs(beta) < smlnum) {
            //
            //           XNORM, BETA may be inaccurate; scale X and recompute them
            //
            bignum = one / smlnum;
        statement_10:
            knt++;
            Rscal(n - 1, bignum, x, incx);
            beta = beta * bignum;
            alpha = alpha * bignum;
            if ((abs(beta) < smlnum) && (knt < 20)) {
                goto statement_10;
            }
            //
            //           New BETA is at most 1, at least SMLNUM
            //
            xnorm = Rnrm2(n - 1, x, incx);
            beta = sign(Rlapy2(alpha, xnorm), alpha);
        }
        savealpha = alpha;
        alpha += beta;
        if (beta < zero) {
            beta = -beta;
            tau = -alpha / beta;
        } else {
            alpha = xnorm * (xnorm / alpha);
            tau = alpha / beta;
            alpha = -alpha;
        }
        //
        if (abs(tau) <= smlnum) {
            //
            //           In the case where the computed TAU ends up being a denormalized number,
            //           it loses relative accuracy. This is a BIG problem. Solution: flush TAU
            //           to ZERO. This explains the next IF statement.
            //
            //           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
            //           (Thanks Pat. Thanks MathWorks.)
            //
            if (savealpha >= zero) {
                tau = zero;
            } else {
                tau = two;
                for (j = 1; j <= n - 1; j = j + 1) {
                    x[(1 + (j - 1) * incx) - 1] = 0.0;
                }
                beta = -savealpha;
            }
            //
        } else {
            //
            //           This is the general case.
            //
            Rscal(n - 1, one / alpha, x, incx);
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
    //     End of Rlarfgp
    //
}

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

INTEGER Rlaneg(INTEGER const n, REAL *d, REAL *lld, REAL const sigma, REAL const /* pivmin */, INTEGER const r) {
    INTEGER return_value = 0;
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
    //     Some architectures propagate Infinities and NaNs very slowly, so
    //     the code computes counts in BLKLEN chunks.  Then a NaN can
    //     propagate at most BLKLEN columns before being detected.  This is
    //     not a general tuning parameter; it needs only to be just large
    //     enough that the overhead is tiny in common cases.
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER negcnt = 0;
    //
    //     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
    REAL t = -sigma;
    INTEGER bj = 0;
    const INTEGER blklen = 128;
    INTEGER neg1 = 0;
    REAL bsav = 0.0;
    INTEGER j = 0;
    REAL dplus = 0.0;
    const REAL zero = 0.0;
    REAL tmp = 0.0;
    bool sawnan = false;
    const REAL one = 1.0;
    for (bj = 1; bj <= r - 1; bj = bj + blklen) {
        neg1 = 0;
        bsav = t;
        for (j = bj; j <= min(bj + blklen - 1, r - 1); j = j + 1) {
            dplus = d[j - 1] + t;
            if (dplus < zero) {
                neg1++;
            }
            tmp = t / dplus;
            t = tmp * lld[j - 1] - sigma;
        }
        sawnan = Risnan(t);
        //     Run a slower version of the above loop if a NaN is detected.
        //     A NaN should occur only with a zero pivot after an infinite
        //     pivot.  In that case, substituting 1 for T/DPLUS is the
        //     correct limit.
        if (sawnan) {
            neg1 = 0;
            t = bsav;
            for (j = bj; j <= min(bj + blklen - 1, r - 1); j = j + 1) {
                dplus = d[j - 1] + t;
                if (dplus < zero) {
                    neg1++;
                }
                tmp = t / dplus;
                if (Risnan(tmp)) {
                    tmp = one;
                }
                t = tmp * lld[j - 1] - sigma;
            }
        }
        negcnt += neg1;
    }
    //
    //     II) lower part: L D L^T - SIGMA I = U- D- U-^T
    REAL p = d[n - 1] - sigma;
    INTEGER neg2 = 0;
    REAL dminus = 0.0;
    for (bj = n - 1; bj >= r; bj = bj - blklen) {
        neg2 = 0;
        bsav = p;
        for (j = bj; j >= max(bj - blklen + 1, r); j = j - 1) {
            dminus = lld[j - 1] + p;
            if (dminus < zero) {
                neg2++;
            }
            tmp = p / dminus;
            p = tmp * d[j - 1] - sigma;
        }
        sawnan = Risnan(p);
        //     As above, run a slower version that substitutes 1 for Inf/Inf.
        //
        if (sawnan) {
            neg2 = 0;
            p = bsav;
            for (j = bj; j >= max(bj - blklen + 1, r); j = j - 1) {
                dminus = lld[j - 1] + p;
                if (dminus < zero) {
                    neg2++;
                }
                tmp = p / dminus;
                if (Risnan(tmp)) {
                    tmp = one;
                }
                p = tmp * d[j - 1] - sigma;
            }
        }
        negcnt += neg2;
    }
    //
    //     III) Twist index
    //       T was shifted by SIGMA initially.
    REAL gamma = (t + sigma) + p;
    if (gamma < zero) {
        negcnt++;
    }
    //
    return_value = negcnt;
    return return_value;
}

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

void Rsb2st_kernels(const char *uplo, bool const &wantz, INTEGER const &ttype, INTEGER const &st, INTEGER const &ed, INTEGER const &sweep, INTEGER const &n, INTEGER const &nb, INTEGER const &ib, REAL *a, INTEGER const &lda, REAL *v, REAL *tau, INTEGER const &ldvt, REAL *work) {
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
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER ajeter = ib + ldvt;
    bool upper = Mlsame(uplo, "U");
    //
    INTEGER dpos = 0;
    INTEGER ofdpos = 0;
    if (upper) {
        dpos = 2 * nb + 1;
        ofdpos = 2 * nb;
    } else {
        dpos = 1;
        ofdpos = 2;
    }
    //
    //     Upper case
    //
    INTEGER vpos = 0;
    INTEGER taupos = 0;
    INTEGER lm = 0;
    const REAL one = 1.0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL ctmp = 0.0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    INTEGER ln = 0;
    if (upper) {
        //
        if (wantz) {
            vpos = mod(sweep - 1, 2) * n + st;
            taupos = mod(sweep - 1, 2) * n + st;
        } else {
            vpos = mod(sweep - 1, 2) * n + st;
            taupos = mod(sweep - 1, 2) * n + st;
        }
        //
        if (ttype == 1) {
            lm = ed - st + 1;
            //
            v[vpos - 1] = one;
            for (i = 1; i <= lm - 1; i = i + 1) {
                v[(vpos + i) - 1] = (a[((ofdpos - i) - 1) + ((st + i) - 1) * lda]);
                a[((ofdpos - i) - 1) + ((st + i) - 1) * lda] = zero;
            }
            ctmp = (a[(ofdpos - 1) + (st - 1) * lda]);
            Rlarfg(lm, ctmp, v[(vpos + 1) - 1], 1, tau[taupos - 1]);
            a[(ofdpos - 1) + (st - 1) * lda] = ctmp;
            //
            lm = ed - st + 1;
            Rlarfy(uplo, lm, v[vpos - 1], 1, (tau[taupos - 1]), a[(dpos - 1) + (st - 1) * lda], lda - 1, work);
        }
        //
        if (ttype == 3) {
            //
            lm = ed - st + 1;
            Rlarfy(uplo, lm, v[vpos - 1], 1, (tau[taupos - 1]), a[(dpos - 1) + (st - 1) * lda], lda - 1, work);
        }
        //
        if (ttype == 2) {
            j1 = ed + 1;
            j2 = min(ed + nb, n);
            ln = ed - st + 1;
            lm = j2 - j1 + 1;
            if (lm > 0) {
                Rlarfx("Left", ln, lm, v[vpos - 1], (tau[taupos - 1]), a[((dpos - nb) - 1) + (j1 - 1) * lda], lda - 1, work);
                //
                if (wantz) {
                    vpos = mod(sweep - 1, 2) * n + j1;
                    taupos = mod(sweep - 1, 2) * n + j1;
                } else {
                    vpos = mod(sweep - 1, 2) * n + j1;
                    taupos = mod(sweep - 1, 2) * n + j1;
                }
                //
                v[vpos - 1] = one;
                for (i = 1; i <= lm - 1; i = i + 1) {
                    v[(vpos + i) - 1] = (a[((dpos - nb - i) - 1) + ((j1 + i) - 1) * lda]);
                    a[((dpos - nb - i) - 1) + ((j1 + i) - 1) * lda] = zero;
                }
                ctmp = (a[((dpos - nb) - 1) + (j1 - 1) * lda]);
                Rlarfg(lm, ctmp, v[(vpos + 1) - 1], 1, tau[taupos - 1]);
                a[((dpos - nb) - 1) + (j1 - 1) * lda] = ctmp;
                //
                Rlarfx("Right", ln - 1, lm, v[vpos - 1], tau[taupos - 1], a[((dpos - nb + 1) - 1) + (j1 - 1) * lda], lda - 1, work);
            }
        }
        //
        //     Lower case
        //
    } else {
        //
        if (wantz) {
            vpos = mod(sweep - 1, 2) * n + st;
            taupos = mod(sweep - 1, 2) * n + st;
        } else {
            vpos = mod(sweep - 1, 2) * n + st;
            taupos = mod(sweep - 1, 2) * n + st;
        }
        //
        if (ttype == 1) {
            lm = ed - st + 1;
            //
            v[vpos - 1] = one;
            for (i = 1; i <= lm - 1; i = i + 1) {
                v[(vpos + i) - 1] = a[((ofdpos + i) - 1) + ((st - 1) - 1) * lda];
                a[((ofdpos + i) - 1) + ((st - 1) - 1) * lda] = zero;
            }
            Rlarfg(lm, a[(ofdpos - 1) + ((st - 1) - 1) * lda], v[(vpos + 1) - 1], 1, tau[taupos - 1]);
            //
            lm = ed - st + 1;
            //
            Rlarfy(uplo, lm, v[vpos - 1], 1, (tau[taupos - 1]), a[(dpos - 1) + (st - 1) * lda], lda - 1, work);
            //
        }
        //
        if (ttype == 3) {
            lm = ed - st + 1;
            //
            Rlarfy(uplo, lm, v[vpos - 1], 1, (tau[taupos - 1]), a[(dpos - 1) + (st - 1) * lda], lda - 1, work);
            //
        }
        //
        if (ttype == 2) {
            j1 = ed + 1;
            j2 = min(ed + nb, n);
            ln = ed - st + 1;
            lm = j2 - j1 + 1;
            //
            if (lm > 0) {
                Rlarfx("Right", lm, ln, v[vpos - 1], tau[taupos - 1], a[((dpos + nb) - 1) + (st - 1) * lda], lda - 1, work);
                //
                if (wantz) {
                    vpos = mod(sweep - 1, 2) * n + j1;
                    taupos = mod(sweep - 1, 2) * n + j1;
                } else {
                    vpos = mod(sweep - 1, 2) * n + j1;
                    taupos = mod(sweep - 1, 2) * n + j1;
                }
                //
                v[vpos - 1] = one;
                for (i = 1; i <= lm - 1; i = i + 1) {
                    v[(vpos + i) - 1] = a[((dpos + nb + i) - 1) + (st - 1) * lda];
                    a[((dpos + nb + i) - 1) + (st - 1) * lda] = zero;
                }
                Rlarfg(lm, a[((dpos + nb) - 1) + (st - 1) * lda], v[(vpos + 1) - 1], 1, tau[taupos - 1]);
                //
                Rlarfx("Left", lm, ln - 1, v[vpos - 1], (tau[taupos - 1]), a[((dpos + nb - 1) - 1) + ((st + 1) - 1) * lda], lda - 1, work);
                //
            }
        }
    }
    //
    //     END OF Rsb2st_kernels
    //
}

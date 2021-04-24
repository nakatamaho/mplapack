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

void Rlagge(INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL *d, REAL *a, INTEGER const lda, INTEGER *iseed, REAL *work, INTEGER &info) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kl < 0 || kl > m - 1) {
        info = -3;
    } else if (ku < 0 || ku > n - 1) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -7;
    }
    if (info < 0) {
        Mxerbla("Rlagge", -info);
        return;
    }
    //
    //     initialize A to diagonal matrix
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
    }
    for (i = 1; i <= min(m, n); i = i + 1) {
        a[(i - 1) + (i - 1) * lda] = d[i - 1];
    }
    //
    //     Quick exit if the user wants a diagonal matrix
    //
    if ((kl == 0) && (ku == 0)) {
        return;
    }
    //
    //     pre- and post-multiply A by random orthogonal matrices
    //
    REAL wn = 0.0;
    REAL wa = 0.0;
    REAL tau = 0.0;
    REAL wb = 0.0;
    const REAL one = 1.0;
    for (i = min(m, n); i >= 1; i = i - 1) {
        if (i < m) {
            //
            //           generate random reflection
            //
            dlarnv(3, iseed, m - i + 1, work);
            wn = Rnrm2(m - i + 1, work, 1);
            wa = sign(wn, &work[1 - 1]);
            if (wn == zero) {
                tau = zero;
            } else {
                wb = work[1 - 1] + wa;
                Rscal(m - i, one / wb, &work[2 - 1], 1);
                work[1 - 1] = one;
                tau = wb / wa;
            }
            //
            //           multiply A(i:m,i:n) by random reflection from the left
            //
            Rgemv("Transpose", m - i + 1, n - i + 1, one, &a[(i - 1) + (i - 1) * lda], lda, work, 1, zero, &work[(m + 1) - 1], 1);
            Rger(m - i + 1, n - i + 1, -tau, work, 1, &work[(m + 1) - 1], 1, &a[(i - 1) + (i - 1) * lda], lda);
        }
        if (i < n) {
            //
            //           generate random reflection
            //
            dlarnv(3, iseed, n - i + 1, work);
            wn = Rnrm2(n - i + 1, work, 1);
            wa = sign(wn, &work[1 - 1]);
            if (wn == zero) {
                tau = zero;
            } else {
                wb = work[1 - 1] + wa;
                Rscal(n - i, one / wb, &work[2 - 1], 1);
                work[1 - 1] = one;
                tau = wb / wa;
            }
            //
            //           multiply A(i:m,i:n) by random reflection from the right
            //
            Rgemv("No transpose", m - i + 1, n - i + 1, one, &a[(i - 1) + (i - 1) * lda], lda, work, 1, zero, &work[(n + 1) - 1], 1);
            Rger(m - i + 1, n - i + 1, -tau, &work[(n + 1) - 1], 1, work, 1, &a[(i - 1) + (i - 1) * lda], lda);
        }
    }
    //
    //     Reduce number of subdiagonals to KL and number of superdiagonals
    //     to KU
    //
    for (i = 1; i <= max(m - 1 - kl, n - 1 - ku); i = i + 1) {
        if (kl <= ku) {
            //
            //           annihilate subdiagonal elements first (necessary if KL = 0)
            //
            if (i <= min(m - 1 - kl, n)) {
                //
                //              generate reflection to annihilate A(kl+i+1:m,i)
                //
                wn = Rnrm2(m - kl - i + 1, &a[((kl + i) - 1) + (i - 1) * lda], 1);
                wa = sign(wn, &a[((kl + i) - 1) + (i - 1) * lda]);
                if (wn == zero) {
                    tau = zero;
                } else {
                    wb = a[((kl + i) - 1) + (i - 1) * lda] + wa;
                    Rscal(m - kl - i, one / wb, &a[((kl + i + 1) - 1) + (i - 1) * lda], 1);
                    a[((kl + i) - 1) + (i - 1) * lda] = one;
                    tau = wb / wa;
                }
                //
                //              apply reflection to A(kl+i:m,i+1:n) from the left
                //
                Rgemv("Transpose", m - kl - i + 1, n - i, one, &a[((kl + i) - 1) + ((i + 1) - 1) * lda], lda, &a[((kl + i) - 1) + (i - 1) * lda], 1, zero, work, 1);
                Rger(m - kl - i + 1, n - i, -tau, &a[((kl + i) - 1) + (i - 1) * lda], 1, work, 1, &a[((kl + i) - 1) + ((i + 1) - 1) * lda], lda);
                a[((kl + i) - 1) + (i - 1) * lda] = -wa;
            }
            //
            if (i <= min(n - 1 - ku, m)) {
                //
                //              generate reflection to annihilate A(i,ku+i+1:n)
                //
                wn = Rnrm2(n - ku - i + 1, &a[(i - 1) + ((ku + i) - 1) * lda], lda);
                wa = sign(wn, &a[(i - 1) + ((ku + i) - 1) * lda]);
                if (wn == zero) {
                    tau = zero;
                } else {
                    wb = a[(i - 1) + ((ku + i) - 1) * lda] + wa;
                    Rscal(n - ku - i, one / wb, &a[(i - 1) + ((ku + i + 1) - 1) * lda], lda);
                    a[(i - 1) + ((ku + i) - 1) * lda] = one;
                    tau = wb / wa;
                }
                //
                //              apply reflection to A(i+1:m,ku+i:n) from the right
                //
                Rgemv("No transpose", m - i, n - ku - i + 1, one, &a[((i + 1) - 1) + ((ku + i) - 1) * lda], lda, &a[(i - 1) + ((ku + i) - 1) * lda], lda, zero, work, 1);
                Rger(m - i, n - ku - i + 1, -tau, work, 1, &a[(i - 1) + ((ku + i) - 1) * lda], lda, &a[((i + 1) - 1) + ((ku + i) - 1) * lda], lda);
                a[(i - 1) + ((ku + i) - 1) * lda] = -wa;
            }
        } else {
            //
            //           annihilate superdiagonal elements first (necessary if
            //           KU = 0)
            //
            if (i <= min(n - 1 - ku, m)) {
                //
                //              generate reflection to annihilate A(i,ku+i+1:n)
                //
                wn = Rnrm2(n - ku - i + 1, &a[(i - 1) + ((ku + i) - 1) * lda], lda);
                wa = sign(wn, &a[(i - 1) + ((ku + i) - 1) * lda]);
                if (wn == zero) {
                    tau = zero;
                } else {
                    wb = a[(i - 1) + ((ku + i) - 1) * lda] + wa;
                    Rscal(n - ku - i, one / wb, &a[(i - 1) + ((ku + i + 1) - 1) * lda], lda);
                    a[(i - 1) + ((ku + i) - 1) * lda] = one;
                    tau = wb / wa;
                }
                //
                //              apply reflection to A(i+1:m,ku+i:n) from the right
                //
                Rgemv("No transpose", m - i, n - ku - i + 1, one, &a[((i + 1) - 1) + ((ku + i) - 1) * lda], lda, &a[(i - 1) + ((ku + i) - 1) * lda], lda, zero, work, 1);
                Rger(m - i, n - ku - i + 1, -tau, work, 1, &a[(i - 1) + ((ku + i) - 1) * lda], lda, &a[((i + 1) - 1) + ((ku + i) - 1) * lda], lda);
                a[(i - 1) + ((ku + i) - 1) * lda] = -wa;
            }
            //
            if (i <= min(m - 1 - kl, n)) {
                //
                //              generate reflection to annihilate A(kl+i+1:m,i)
                //
                wn = Rnrm2(m - kl - i + 1, &a[((kl + i) - 1) + (i - 1) * lda], 1);
                wa = sign(wn, &a[((kl + i) - 1) + (i - 1) * lda]);
                if (wn == zero) {
                    tau = zero;
                } else {
                    wb = a[((kl + i) - 1) + (i - 1) * lda] + wa;
                    Rscal(m - kl - i, one / wb, &a[((kl + i + 1) - 1) + (i - 1) * lda], 1);
                    a[((kl + i) - 1) + (i - 1) * lda] = one;
                    tau = wb / wa;
                }
                //
                //              apply reflection to A(kl+i:m,i+1:n) from the left
                //
                Rgemv("Transpose", m - kl - i + 1, n - i, one, &a[((kl + i) - 1) + ((i + 1) - 1) * lda], lda, &a[((kl + i) - 1) + (i - 1) * lda], 1, zero, work, 1);
                Rger(m - kl - i + 1, n - i, -tau, &a[((kl + i) - 1) + (i - 1) * lda], 1, work, 1, &a[((kl + i) - 1) + ((i + 1) - 1) * lda], lda);
                a[((kl + i) - 1) + (i - 1) * lda] = -wa;
            }
        }
        //
        if (i <= n) {
            for (j = kl + i + 1; j <= m; j = j + 1) {
                a[(j - 1) + (i - 1) * lda] = zero;
            }
        }
        //
        if (i <= m) {
            for (j = ku + i + 1; j <= n; j = j + 1) {
                a[(i - 1) + (j - 1) * lda] = zero;
            }
        }
    }
    //
    //     End of Rlagge
    //
}

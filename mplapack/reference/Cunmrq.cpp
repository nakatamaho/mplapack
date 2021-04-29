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

void Cunmrq(const char *side, const char *trans, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *c, INTEGER const ldc, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    bool left = Mlsame(side, "L");
    bool notran = Mlsame(trans, "N");
    bool lquery = (lwork == -1);
    //
    //
    INTEGER nq = 0;
    INTEGER nw = 0;
    if (left) {
        nq = m;
        nw = max((INTEGER)1, n);
    } else {
        nq = n;
        nw = max((INTEGER)1, m);
    }
    if (!left && !Mlsame(side, "R")) {
        info = -1;
    } else if (!notran && !Mlsame(trans, "C")) {
        info = -2;
    } else if (m < 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (k < 0 || k > nq) {
        info = -5;
    } else if (lda < max((INTEGER)1, k)) {
        info = -7;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -10;
    } else if (lwork < nw && !lquery) {
        info = -12;
    }
    //
    INTEGER lwkopt = 0;
    const INTEGER nbmax = 64;
    INTEGER nb = 0;
    const INTEGER ldt = nbmax + 1;
    const INTEGER tsize = ldt * nbmax;
    char side_trans[3];
    side_trans[0] = side[0];
    side_trans[1] = trans[0];
    side_trans[2] = '\0';
    if (info == 0) {
        //
        //        Compute the workspace requirements
        //
        if (m == 0 || n == 0) {
            lwkopt = 1;
        } else {
            nb = min({nbmax, iMlaenv(1, "Cunmrq", side_trans, m, n, k, -1)});
            lwkopt = nw * nb + tsize;
        }
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Cunmrq", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    INTEGER nbmin = 2;
    INTEGER ldwork = nw;
    if (nb > 1 && nb < k) {
        if (lwork < nw * nb + tsize) {
            nb = (lwork - tsize) / ldwork;
            nbmin = max({(INTEGER)2, iMlaenv(2, "Cunmrq", side_trans, m, n, k, -1)});
        }
    }
    //
    INTEGER iinfo = 0;
    INTEGER iwt = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    INTEGER ni = 0;
    INTEGER mi = 0;
    char transt;
    INTEGER i = 0;
    INTEGER ib = 0;
    if (nb < nbmin || nb >= k) {
        //
        //        Use unblocked code
        //
        Cunmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, iinfo);
    } else {
        //
        //        Use blocked code
        //
        iwt = 1 + nw * nb;
        if ((left && !notran) || (!left && notran)) {
            i1 = 1;
            i2 = k;
            i3 = nb;
        } else {
            i1 = ((k - 1) / nb) * nb + 1;
            i2 = 1;
            i3 = -nb;
        }
        //
        if (left) {
            ni = n;
        } else {
            mi = m;
        }
        //
        if (notran) {
            transt = 'C';
        } else {
            transt = 'N';
        }
        //
        for (i = i1; i <= i2; i = i + i3) {
            ib = min(nb, k - i + 1);
            //
            //           Form the triangular factor of the block reflector
            //           H = H(i+ib-1) . . . H(i+1) H(i)
            //
            Clarft("Backward", "Rowwise", nq - k + i + ib - 1, ib, &a[(i - 1)], lda, &tau[i - 1], &work[iwt - 1], ldt);
            if (left) {
                //
                //              H or H**H is applied to C(1:m-k+i+ib-1,1:n)
                //
                mi = m - k + i + ib - 1;
            } else {
                //
                //              H or H**H is applied to C(1:m,1:n-k+i+ib-1)
                //
                ni = n - k + i + ib - 1;
            }
            //
            //           Apply H or H**H
            //
            Clarfb(side, &transt, "Backward", "Rowwise", mi, ni, ib, &a[(i - 1)], lda, &work[iwt - 1], ldt, c, ldc, work, ldwork);
        }
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Cunmrq
    //
}

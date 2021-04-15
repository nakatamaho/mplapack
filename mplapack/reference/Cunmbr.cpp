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

void Cunmbr(const char *vect, const char *side, const char *trans, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *c, INTEGER const ldc, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    bool applyq = Mlsame(vect, "Q");
    bool left = Mlsame(side, "L");
    bool notran = Mlsame(trans, "N");
    bool lquery = (lwork == -1);
    //
    //
    INTEGER nq = 0;
    INTEGER nw = 0;
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if (m == 0 || n == 0) {
        nw = 0;
    }
    if (!applyq && !Mlsame(vect, "P")) {
        info = -1;
    } else if (!left && !Mlsame(side, "R")) {
        info = -2;
    } else if (!notran && !Mlsame(trans, "C")) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (k < 0) {
        info = -6;
    } else if ((applyq && lda < max((INTEGER)1, nq)) || (!applyq && lda < max((INTEGER)1, min(nq, k)))) {
        info = -8;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -11;
    } else if (lwork < max((INTEGER)1, nw) && !lquery) {
        info = -13;
    }
    //
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    char side_trans[3];
    side_trans[0] = side[0];
    side_trans[1] = trans[0];
    side_trans[2] = '\0';
    if (info == 0) {
        if (nw > 0) {
            if (applyq) {
                if (left) {
                    nb = iMlaenv(1, "Cunmqr", side_trans, m - 1, n, m - 1, -1);
                } else {
                    nb = iMlaenv(1, "Cunmqr", side_trans, m, n - 1, n - 1, -1);
                }
            } else {
                if (left) {
                    nb = iMlaenv(1, "Cunmlq", side_trans, m - 1, n, m - 1, -1);
                } else {
                    nb = iMlaenv(1, "Cunmlq", side_trans, m, n - 1, n - 1, -1);
                }
            }
            lwkopt = max((INTEGER)1, nw * nb);
        } else {
            lwkopt = 1;
        }
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Cunmbr", -info);
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
    INTEGER iinfo = 0;
    INTEGER mi = 0;
    INTEGER ni = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    char transt;
    if (applyq) {
        //
        //        Apply Q
        //
        if (nq >= k) {
            //
            //           Q was determined by a call to Cgebrd with nq >= k
            //
            Cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, iinfo);
        } else if (nq > 1) {
            //
            //           Q was determined by a call to Cgebrd with nq < k
            //
            if (left) {
                mi = m - 1;
                ni = n;
                i1 = 2;
                i2 = 1;
            } else {
                mi = m;
                ni = n - 1;
                i1 = 1;
                i2 = 2;
            }
            Cunmqr(side, trans, mi, ni, nq - 1, &a[(2 - 1)], lda, tau, &c[(i1 - 1) + (i2 - 1) * ldc], ldc, work, lwork, iinfo);
        }
    } else {
        //
        //        Apply P
        //
        if (notran) {
            transt = 'C';
        } else {
            transt = 'N';
        }
        if (nq > k) {
            //
            //           P was determined by a call to Cgebrd with nq > k
            //
            Cunmlq(side, &transt, m, n, k, a, lda, tau, c, ldc, work, lwork, iinfo);
        } else if (nq > 1) {
            //
            //           P was determined by a call to Cgebrd with nq <= k
            //
            if (left) {
                mi = m - 1;
                ni = n;
                i1 = 2;
                i2 = 1;
            } else {
                mi = m;
                ni = n - 1;
                i1 = 1;
                i2 = 2;
            }
            Cunmlq(side, &transt, mi, ni, nq - 1, &a[(2 - 1) * lda], lda, tau, &c[(i1 - 1) + (i2 - 1) * ldc], ldc, work, lwork, iinfo);
        }
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Cunmbr
    //
}

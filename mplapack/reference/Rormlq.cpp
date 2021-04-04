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

void Rormlq(const char *side, const char *trans, INTEGER const &m, INTEGER const &n, INTEGER const &k, REAL *a, INTEGER const &lda, REAL *tau, REAL *c, INTEGER const &ldc, REAL *work, INTEGER const &lwork, INTEGER &info) {
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
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if (!left && !Mlsame(side, "R")) {
        info = -1;
    } else if (!notran && !Mlsame(trans, "T")) {
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
    } else if (lwork < max((INTEGER)1, nw) && !lquery) {
        info = -12;
    }
    //
    const INTEGER nbmax = 64;
    INTEGER nb = 0;
    const INTEGER ldt = nbmax + 1;
    const INTEGER tsize = ldt * nbmax;
    INTEGER lwkopt = 0;
    if (info == 0) {
        //
        //        Compute the workspace requirements
        //
        nb = min(nbmax, iMlaenv[("Rormlq" - 1) * ldiMlaenv]);
        lwkopt = max((INTEGER)1, nw) * nb + tsize;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rormlq", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0 || k == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    INTEGER nbmin = 2;
    INTEGER ldwork = nw;
    if (nb > 1 && nb < k) {
        if (lwork < nw * nb + tsize) {
            nb = (lwork - tsize) / ldwork;
            nbmin = max(2, iMlaenv[(2 - 1) + ("Rormlq" - 1) * ldiMlaenv]);
        }
    }
    //
    INTEGER iinfo = 0;
    INTEGER iwt = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    INTEGER ni = 0;
    INTEGER jc = 0;
    INTEGER mi = 0;
    INTEGER ic = 0;
    str<1> transt = char0;
    INTEGER i = 0;
    INTEGER ib = 0;
    if (nb < nbmin || nb >= k) {
        //
        //        Use unblocked code
        //
        Rorml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, iinfo);
    } else {
        //
        //        Use blocked code
        //
        iwt = 1 + nw * nb;
        if ((left && notran) || (!left && !notran)) {
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
            jc = 1;
        } else {
            mi = m;
            ic = 1;
        }
        //
        if (notran) {
            transt = "T";
        } else {
            transt = "N";
        }
        //
        for (i = i1; i <= i2; i = i + i3) {
            ib = min(nb, k - i + 1);
            //
            //           Form the triangular factor of the block reflector
            //           H = H(i) H(i+1) . . . H(i+ib-1)
            //
            Rlarft("Forward", "Rowwise", nq - i + 1, ib, a[(i - 1) + (i - 1) * lda], lda, tau[i - 1], work[iwt - 1], ldt);
            if (left) {
                //
                //              H or H**T is applied to C(i:m,1:n)
                //
                mi = m - i + 1;
                ic = i;
            } else {
                //
                //              H or H**T is applied to C(1:m,i:n)
                //
                ni = n - i + 1;
                jc = i;
            }
            //
            //           Apply H or H**T
            //
            Rlarfb(side, transt, "Forward", "Rowwise", mi, ni, ib, a[(i - 1) + (i - 1) * lda], lda, work[iwt - 1], ldt, c[(ic - 1) + (jc - 1) * ldc], ldc, work, ldwork);
        }
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Rormlq
    //
}

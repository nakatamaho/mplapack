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

void Rgebd2(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *d, REAL *e, REAL *tauq, REAL *taup, REAL *work, INTEGER &info) {
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
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    if (info < 0) {
        Mxerbla("Rgebd2", -info);
        return;
    }
    //
    INTEGER i = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    if (m >= n) {
        //
        //        Reduce to upper bidiagonal form
        //
        for (i = 1; i <= n; i = i + 1) {
            //
            //           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
            //
	    Rlarfg(m - i + 1, a[(i - 1) + (i - 1) * lda], &a[(min(i+1,m)-1) + (i-1) * lda], 1, tauq[i - 1]);
            d[i - 1] = a[(i - 1) + (i - 1) * lda];
            a[(i - 1) + (i - 1) * lda] = one;
            //
            //           Apply H(i) to A(i:m,i+1:n) from the left
            //
            if (i < n) {
                Rlarf("Left", m - i + 1, n - i, &a[(i - 1) + (i - 1) * lda], 1, tauq[i - 1], &a[(i - 1) + ((i + 1) - 1) * lda], lda, work);
            }
            a[(i - 1) + (i - 1) * lda] = d[i - 1];
            //
            if (i < n) {
                //
                //              Generate elementary reflector G(i) to annihilate
                //              A(i,i+2:n)
                //
	      Rlarfg(n - i, a[(i-1)+((i + 1)-1)*lda], &a[(i-1)+(min(i + 2,n)-1)*lda], lda, taup[i-1]);
        e[i-1] = a[(i-1)+((i + 1)-1)*lda];
        a[(i-1)+((i + 1)-1)*lda] = one;
        //
        //              Apply G(i) to A(i+1:m,i+1:n) from the right
        //
        Rlarf("Right", m - i, n - i, &a[(i-1)+((i + 1)-1)*lda], lda, taup[i-1], &a[((i + 1)-1)+((i + 1)-1)*lda], lda, work);
        a[(i-1)+((i + 1)-1)*lda] = e[i-1];
            } else {
                taup[i - 1] = zero;
            }
        }
    } else {
        //
        //        Reduce to lower bidiagonal form
        //
        for (i = 1; i <= m; i = i + 1) {
            //
            //           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
            //
	  Rlarfg(n - i + 1, a[(i-1)+(i-1)*lda], &a[(i-1)+(min(i + 1, n)-1)*lda], lda, taup[i-1]);
      d[i-1] = a[(i-1)+(i-1)*lda];
      a[(i-1)+(i-1)*lda] = one;
      //
      //           Apply G(i) to A(i+1:m,i:n) from the right
      //
      if (i < m) {
                Rlarf("Right", m - i, n - i + 1, &a[(i - 1) + (i - 1) * lda], lda, taup[i - 1], &a[((i + 1) - 1) + (i - 1) * lda], lda, work);
      }
      a[(i-1)+(i-1)*lda] = d[i-1];
      //
      if (i < m) {
                //
                //              Generate elementary reflector H(i) to annihilate
                //              A(i+2:m,i)
                //
	Rlarfg(m - i, a[((i + 1) - 1) + (i - 1) * lda], &a[(min(i + 2, m)-1) + (i-1)* lda], 1, tauq[i - 1]);
                e[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                a[((i + 1) - 1) + (i - 1) * lda] = one;
                //
                //              Apply H(i) to A(i+1:m,i+1:n) from the left
                //
                Rlarf("Left", m - i, n - i, &a[((i + 1) - 1) + (i - 1) * lda], 1, tauq[i - 1], &a[((i + 1) - 1) + ((i + 1) - 1) * lda], lda, work);
                a[((i + 1) - 1) + (i - 1) * lda] = e[i - 1];
      }
      else {
                tauq[i - 1] = zero;
      }
        }
    }
    //
    //     End of Rgebd2
    //
}

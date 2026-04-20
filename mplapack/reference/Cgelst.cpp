/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZGELST.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Cgelst(const char *trans, INTEGER const m, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
    //
    // Test the input arguments.
    //
    info = 0;
    INTEGER mn = min(m, n);
    bool lquery = (lwork == -1);
    if (!(Mlsame(trans, "N") || Mlsame(trans, "C"))) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (nrhs < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -6;
    } else if (ldb < max((INTEGER)1, m, n)) {
        info = -8;
    } else if (lwork < max((INTEGER)1, mn + max(mn, nrhs)) && !lquery) {
        info = -10;
    }
    //
    // Figure out optimal block size and optimal workspace size
    //
    bool tpsd = false;
    INTEGER nb = 0;
    INTEGER mnnrhs = 0;
    INTEGER lwopt = 0;
    if (info == 0 || info == -10) {
        //
        tpsd = true;
        if (Mlsame(trans, "N")) {
            tpsd = false;
        }
        //
        nb = iMlaenv(1, "Cgelst", " ", m, n, -1, -1);
        //
        mnnrhs = max(mn, nrhs);
        lwopt = max((INTEGER)1, (mn + mnnrhs) * nb);
        work[1 - 1] = castREAL(lwopt);
        //
    }
    //
    if (info != 0) {
        Mxerbla("Cgelst", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    // Quick return if possible
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (min(m, n, nrhs) == 0) {
        Claset("Full", max(m, n), nrhs, czero, czero, b, ldb);
        work[1 - 1] = castREAL(lwopt);
        return;
    }
    //
    // *GEQRT and *GELQT routines cannot accept NB larger than min(M,N)
    //
    if (nb > mn) {
        nb = mn;
    }
    //
    // Determine the block size from the supplied LWORK
    // ( at this stage we know that LWORK >= (minimum required workspace,
    // but it may be less than optimal)
    //
    nb = min(nb, lwork / (mn + mnnrhs));
    //
    // The minimum value of NB, when blocked code is used
    //
    INTEGER nbmin = max((INTEGER)2, iMlaenv(2, "Cgelst", " ", m, n, -1, -1));
    //
    if (nb < nbmin) {
        nb = 1;
    }
    //
    // Get machine parameters
    //
    REAL smlnum = Rlamch("S") / Rlamch("P");
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    // Scale A, B if max element outside range [SMLNUM,BIGNUM]
    //
    REAL rwork[1];
    REAL anrm = Clange("M", m, n, a, lda, rwork);
    INTEGER iascl = 0;
    const REAL zero = 0.0;
    if (anrm > zero && anrm < smlnum) {
        //
        // Scale matrix norm up to SMLNUM
        //
        Clascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        iascl = 1;
    } else if (anrm > bignum) {
        //
        // Scale matrix norm down to BIGNUM
        //
        Clascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
        iascl = 2;
    } else if (anrm == zero) {
        //
        // Matrix all zero. Return zero solution.
        //
        Claset("Full", max(m, n), nrhs, czero, czero, b, ldb);
        work[1 - 1] = castREAL(lwopt);
        return;
    }
    //
    INTEGER brow = m;
    if (tpsd) {
        brow = n;
    }
    REAL bnrm = Clange("M", brow, nrhs, b, ldb, rwork);
    INTEGER ibscl = 0;
    if (bnrm > zero && bnrm < smlnum) {
        //
        // Scale matrix norm up to SMLNUM
        //
        Clascl("G", 0, 0, bnrm, smlnum, brow, nrhs, b, ldb, info);
        ibscl = 1;
    } else if (bnrm > bignum) {
        //
        // Scale matrix norm down to BIGNUM
        //
        Clascl("G", 0, 0, bnrm, bignum, brow, nrhs, b, ldb, info);
        ibscl = 2;
    }
    //
    INTEGER scllen = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    if (m >= n) {
        //
        // M > N:
        // Compute the blocked QR factorization of A,
        // using the compact WY representation of Q,
        // workspace at least N, optimally N*NB.
        //
        Cgeqrt(m, n, nb, a, lda, &work[1 - 1], nb, &work[(mn * nb + 1) - 1], info);
        //
        if (!tpsd) {
            //
            // M > N, A is not transposed:
            // Overdetermined system of equations,
            // least-squares problem, min || A * X - B ||.
            //
            // Compute B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.
            //
            Cgemqrt("Left", "Conjugate transpose", m, nrhs, n, nb, a, lda, &work[1 - 1], nb, b, ldb, &work[(mn * nb + 1) - 1], info);
            //
            // Compute B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
            //
            Ctrtrs("Upper", "No transpose", "Non-unit", n, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            scllen = n;
            //
        } else {
            //
            // M > N, A is transposed:
            // Underdetermined system of equations,
            // minimum norm solution of A**T * X = B.
            //
            // Compute B := inv(R**T) * B in two row blocks of B.
            //
            // Block 1: B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
            //
            Ctrtrs("Upper", "Conjugate transpose", "Non-unit", n, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            // Block 2: Zero out all rows below the N-th row in B:
            // B(N+1:M,1:NRHS) = ZERO
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                for (i = n + 1; i <= m; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = zero;
                }
            }
            //
            // Compute B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.
            //
            Cgemqrt("Left", "No transpose", m, nrhs, n, nb, a, lda, &work[1 - 1], nb, b, ldb, &work[(mn * nb + 1) - 1], info);
            //
            scllen = m;
            //
        }
        //
    } else {
        //
        // M < N:
        // Compute the blocked LQ factorization of A,
        // using the compact WY representation of Q,
        // workspace at least M, optimally M*NB.
        //
        Cgelqt(m, n, nb, a, lda, &work[1 - 1], nb, &work[(mn * nb + 1) - 1], info);
        //
        if (!tpsd) {
            //
            // M < N, A is not transposed:
            // Underdetermined system of equations,
            // minimum norm solution of A * X = B.
            //
            // Compute B := inv(L) * B in two row blocks of B.
            //
            // Block 1: B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
            //
            Ctrtrs("Lower", "No transpose", "Non-unit", m, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            // Block 2: Zero out all rows below the M-th row in B:
            // B(M+1:N,1:NRHS) = ZERO
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                for (i = m + 1; i <= n; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = zero;
                }
            }
            //
            // Compute B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.
            //
            Cgemlqt("Left", "Conjugate transpose", n, nrhs, m, nb, a, lda, &work[1 - 1], nb, b, ldb, &work[(mn * nb + 1) - 1], info);
            //
            scllen = n;
            //
        } else {
            //
            // M < N, A is transposed:
            // Overdetermined system of equations,
            // least-squares problem, min || A**T * X - B ||.
            //
            // Compute B(1:N,1:NRHS) := Q * B(1:N,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.
            //
            Cgemlqt("Left", "No transpose", n, nrhs, m, nb, a, lda, &work[1 - 1], nb, b, ldb, &work[(mn * nb + 1) - 1], info);
            //
            // Compute B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
            //
            Ctrtrs("Lower", "Conjugate transpose", "Non-unit", m, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            scllen = m;
            //
        }
        //
    }
    //
    // Undo scaling
    //
    if (iascl == 1) {
        Clascl("G", 0, 0, anrm, smlnum, scllen, nrhs, b, ldb, info);
    } else if (iascl == 2) {
        Clascl("G", 0, 0, anrm, bignum, scllen, nrhs, b, ldb, info);
    }
    if (ibscl == 1) {
        Clascl("G", 0, 0, smlnum, bnrm, scllen, nrhs, b, ldb, info);
    } else if (ibscl == 2) {
        Clascl("G", 0, 0, bignum, bnrm, scllen, nrhs, b, ldb, info);
    }
    //
    work[1 - 1] = castREAL(lwopt);
    //
    // End of Cgelst
    //
}

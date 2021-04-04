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

void Rsbev_2stage(const char *jobz, const char *uplo, INTEGER const &n, INTEGER const &kd, REAL *ab, INTEGER const &ldab, REAL *w, REAL *z, INTEGER const &ldz, REAL *work, INTEGER const &lwork, INTEGER &info) {
    //
    //  -- LAPACK driver routine --
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
    //     Test the input parameters.
    //
    bool wantz = Mlsame(jobz, "V");
    bool lower = Mlsame(uplo, "L");
    bool lquery = (lwork == -1);
    //
    info = 0;
    if (!(Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(lower || Mlsame(uplo, "U"))) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (kd < 0) {
        info = -4;
    } else if (ldab < kd + 1) {
        info = -6;
    } else if (ldz < 1 || (wantz && ldz < n)) {
        info = -9;
    }
    //
    INTEGER lwmin = 0;
    INTEGER ib = 0;
    INTEGER lhtrd = 0;
    INTEGER lwtrd = 0;
    if (info == 0) {
        if (n <= 1) {
            lwmin = 1;
            work[1 - 1] = lwmin;
        } else {
            ib = iMlaenv2stage[(2 - 1) + ("Rsytrd_SB2ST" - 1) * ldiMlaenv2stage];
            lhtrd = iMlaenv2stage[(3 - 1) + ("Rsytrd_SB2ST" - 1) * ldiMlaenv2stage];
            lwtrd = iMlaenv2stage[(4 - 1) + ("Rsytrd_SB2ST" - 1) * ldiMlaenv2stage];
            lwmin = n + lhtrd + lwtrd;
            work[1 - 1] = lwmin;
        }
        //
        if (lwork < lwmin && !lquery) {
            info = -11;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rsbev_2stage ", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    const REAL one = 1.0;
    if (n == 1) {
        if (lower) {
            w[1 - 1] = ab[(1 - 1)];
        } else {
            w[1 - 1] = ab[((kd + 1) - 1)];
        }
        if (wantz) {
            z[(1 - 1)] = one;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    REAL safmin = dlamch("Safe minimum");
    REAL eps = dlamch("Precision");
    REAL smlnum = safmin / eps;
    REAL bignum = one / smlnum;
    REAL rmin = sqrt(smlnum);
    REAL rmax = sqrt(bignum);
    //
    //     Scale matrix to allowable range, if necessary.
    //
    REAL anrm = Rlansb[("M" - 1) + (uplo - 1) * ldRlansb];
    INTEGER iscale = 0;
    const REAL zero = 0.0;
    REAL sigma = 0.0;
    if (anrm > zero && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        if (lower) {
            Rlascl("B", kd, kd, one, sigma, n, n, ab, ldab, info);
        } else {
            Rlascl("Q", kd, kd, one, sigma, n, n, ab, ldab, info);
        }
    }
    //
    //     Call Rsytrd_SB2ST to reduce symmetric band matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indhous = inde + n;
    INTEGER indwrk = indhous + lhtrd;
    INTEGER llwork = lwork - indwrk + 1;
    //
    INTEGER iinfo = 0;
    Rsytrd_sb2st("N", jobz, uplo, n, kd, ab, ldab, w, work[inde - 1], work[indhous - 1], lhtrd, work[indwrk - 1], llwork, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, call SSTEQR.
    //
    if (!wantz) {
        Rsterf(n, w, work[inde - 1], info);
    } else {
        Rsteqr(jobz, n, w, work[inde - 1], z, ldz, work[indwrk - 1], info);
    }
    //
    //     If matrix was scaled, then rescale eigenvalues appropriately.
    //
    INTEGER imax = 0;
    if (iscale == 1) {
        if (info == 0) {
            imax = n;
        } else {
            imax = info - 1;
        }
        Rscal(imax, one / sigma, w, 1);
    }
    //
    //     Set WORK(1) to optimal workspace size.
    //
    work[1 - 1] = lwmin;
    //
    //     End of Rsbev_2stage
    //
}

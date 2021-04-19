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

void Rsbevd(const char *jobz, const char *uplo, INTEGER const n, INTEGER const kd, REAL *ab, INTEGER const ldab, REAL *w, REAL *z, INTEGER const ldz, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
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
    bool lquery = (lwork == -1 || liwork == -1);
    //
    info = 0;
    INTEGER liwmin = 0;
    INTEGER lwmin = 0;
    if (n <= 1) {
        liwmin = 1;
        lwmin = 1;
    } else {
        if (wantz) {
            liwmin = 3 + 5 * n;
            lwmin = 1 + 5 * n + 2 * pow2(n);
        } else {
            liwmin = 1;
            lwmin = 2 * n;
        }
    }
    if (!(wantz || Mlsame(jobz, "N"))) {
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
    if (info == 0) {
        work[1 - 1] = lwmin;
        iwork[1 - 1] = liwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -11;
        } else if (liwork < liwmin && !lquery) {
            info = -13;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rsbevd", -info);
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
        w[1 - 1] = ab[(1 - 1)];
        if (wantz) {
            z[(1 - 1)] = one;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    REAL safmin = Rlamch("Safe minimum");
    REAL eps = Rlamch("Precision");
    REAL smlnum = safmin / eps;
    REAL bignum = one / smlnum;
    REAL rmin = sqrt(smlnum);
    REAL rmax = sqrt(bignum);
    //
    //     Scale matrix to allowable range, if necessary.
    //
    REAL anrm = Rlansb("M", uplo, n, kd, ab, ldab, work);
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
    //     Call Rsbtrd to reduce symmetric band matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indwrk = inde + n;
    INTEGER indwk2 = indwrk + n * n;
    INTEGER llwrk2 = lwork - indwk2 + 1;
    INTEGER iinfo = 0;
    Rsbtrd(jobz, uplo, n, kd, ab, ldab, w, &work[inde - 1], z, ldz, &work[indwrk - 1], iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, call SSTEDC.
    //
    if (!wantz) {
        Rsterf(n, w, &work[inde - 1], info);
    } else {
        Rstedc("I", n, w, &work[inde - 1], &work[indwrk - 1], n, &work[indwk2 - 1], llwrk2, iwork, liwork, info);
        Rgemm("N", "N", n, n, n, one, z, ldz, &work[indwrk - 1], n, zero, &work[indwk2 - 1], n);
        Rlacpy("A", n, n, &work[indwk2 - 1], n, z, ldz);
    }
    //
    //     If matrix was scaled, then rescale eigenvalues appropriately.
    //
    if (iscale == 1) {
        Rscal(n, one / sigma, w, 1);
    }
    //
    work[1 - 1] = lwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Rsbevd
    //
}

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

void Chbgvx(const char *jobz, const char *range, const char *uplo, INTEGER const &n, INTEGER const &ka, INTEGER const &kb, COMPLEX *ab, INTEGER const &ldab, COMPLEX *bb, INTEGER const &ldbb, COMPLEX *q, INTEGER const &ldq, REAL const &vl, REAL const &vu, INTEGER const &il, INTEGER const &iu, REAL const &abstol, INTEGER &m, REAL *w, COMPLEX *z, INTEGER const &ldz, COMPLEX *work, REAL *rwork, arr_ref<INTEGER> iwork, arr_ref<INTEGER> ifail, INTEGER &info) {
    bool wantz = false;
    bool upper = false;
    bool alleig = false;
    bool valeig = false;
    bool indeig = false;
    INTEGER iinfo = 0;
    INTEGER indd = 0;
    INTEGER inde = 0;
    INTEGER indrwk = 0;
    INTEGER indwrk = 0;
    str<1> vect = char0;
    bool test = false;
    const REAL zero = 0.0;
    INTEGER indee = 0;
    INTEGER i = 0;
    str<1> order = char0;
    INTEGER indibl = 0;
    INTEGER indisp = 0;
    INTEGER indiwk = 0;
    INTEGER nsplit = 0;
    INTEGER j = 0;
    const COMPLEX cone = (1.0, 0.0);
    const COMPLEX czero = (0.0, 0.0);
    REAL tmp1 = 0.0;
    INTEGER jj = 0;
    INTEGER itmp1 = 0;
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
    wantz = Mlsame(jobz, "V");
    upper = Mlsame(uplo, "U");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    //
    info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(alleig || valeig || indeig)) {
        info = -2;
    } else if (!(upper || Mlsame(uplo, "L"))) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (ka < 0) {
        info = -5;
    } else if (kb < 0 || kb > ka) {
        info = -6;
    } else if (ldab < ka + 1) {
        info = -8;
    } else if (ldbb < kb + 1) {
        info = -10;
    } else if (ldq < 1 || (wantz && ldq < n)) {
        info = -12;
    } else {
        if (valeig) {
            if (n > 0 && vu <= vl) {
                info = -14;
            }
        } else if (indeig) {
            if (il < 1 || il > max((INTEGER)1, n)) {
                info = -15;
            } else if (iu < min(n, il) || iu > n) {
                info = -16;
            }
        }
    }
    if (info == 0) {
        if (ldz < 1 || (wantz && ldz < n)) {
            info = -21;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Chbgvx", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    m = 0;
    if (n == 0) {
        return;
    }
    //
    //     Form a split Cholesky factorization of B.
    //
    Cpbstf(uplo, n, kb, bb, ldbb, info);
    if (info != 0) {
        info += n;
        return;
    }
    //
    //     Transform problem to standard eigenvalue problem.
    //
    Chbgst(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, work, rwork, iinfo);
    //
    //     Solve the standard eigenvalue problem.
    //     Reduce Hermitian band matrix to tridiagonal form.
    //
    indd = 1;
    inde = indd + n;
    indrwk = inde + n;
    indwrk = 1;
    if (wantz) {
        vect = "U";
    } else {
        vect = "N";
    }
    Chbtrd(vect, uplo, n, ka, ab, ldab, rwork[indd - 1], rwork[inde - 1], q, ldq, work[indwrk - 1], iinfo);
    //
    //     If all eigenvalues are desired and ABSTOL is less than or equal
    //     to zero, then call Rsterf or Csteqr.  If this fails for some
    //     eigenvalue, then try Rstebz.
    //
    test = false;
    if (indeig) {
        if (il == 1 && iu == n) {
            test = true;
        }
    }
    if ((alleig || test) && (abstol <= zero)) {
        Rcopy(n, rwork[indd - 1], 1, w, 1);
        indee = indrwk + 2 * n;
        Rcopy(n - 1, rwork[inde - 1], 1, rwork[indee - 1], 1);
        if (!wantz) {
            Rsterf(n, w, rwork[indee - 1], info);
        } else {
            Clacpy("A", n, n, q, ldq, z, ldz);
            Csteqr(jobz, n, w, rwork[indee - 1], z, ldz, rwork[indrwk - 1], info);
            if (info == 0) {
                for (i = 1; i <= n; i = i + 1) {
                    ifail[i - 1] = 0;
                }
            }
        }
        if (info == 0) {
            m = n;
            goto statement_30;
        }
        info = 0;
    }
    //
    //     Otherwise, call Rstebz and, if eigenvectors are desired,
    //     call Cstein.
    //
    if (wantz) {
        order = "B";
    } else {
        order = "E";
    }
    indibl = 1;
    indisp = indibl + n;
    indiwk = indisp + n;
    Rstebz(range, order, n, vl, vu, il, iu, abstol, rwork[indd - 1], rwork[inde - 1], m, nsplit, w, iwork[indibl - 1], iwork[indisp - 1], rwork[indrwk - 1], iwork[indiwk - 1], info);
    //
    if (wantz) {
        Cstein(n, rwork[indd - 1], rwork[inde - 1], m, w, iwork[indibl - 1], iwork[indisp - 1], z, ldz, rwork[indrwk - 1], iwork[indiwk - 1], ifail, info);
        //
        //        Apply unitary matrix used in reduction to tridiagonal
        //        form to eigenvectors returned by Cstein.
        //
        for (j = 1; j <= m; j = j + 1) {
            Ccopy(n, z[(j - 1) * ldz], 1, work[1 - 1], 1);
            Cgemv("N", n, n, cone, q, ldq, work, 1, czero, z[(j - 1) * ldz], 1);
        }
    }
//
statement_30:
    //
    //     If eigenvalues are not in order, then sort them, along with
    //     eigenvectors.
    //
    if (wantz) {
        for (j = 1; j <= m - 1; j = j + 1) {
            i = 0;
            tmp1 = w[j - 1];
            for (jj = j + 1; jj <= m; jj = jj + 1) {
                if (w[jj - 1] < tmp1) {
                    i = jj;
                    tmp1 = w[jj - 1];
                }
            }
            //
            if (i != 0) {
                itmp1 = iwork[(indibl + i - 1) - 1];
                w[i - 1] = w[j - 1];
                iwork[(indibl + i - 1) - 1] = iwork[(indibl + j - 1) - 1];
                w[j - 1] = tmp1;
                iwork[(indibl + j - 1) - 1] = itmp1;
                Cswap(n, z[(i - 1) * ldz], 1, z[(j - 1) * ldz], 1);
                if (info != 0) {
                    itmp1 = ifail[i - 1];
                    ifail[i - 1] = ifail[j - 1];
                    ifail[j - 1] = itmp1;
                }
            }
        }
    }
    //
    //     End of Chbgvx
    //
}

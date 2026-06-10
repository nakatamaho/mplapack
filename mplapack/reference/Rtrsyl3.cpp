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

// Derived from LAPACK routine DTRSYL3.
// Original LAPACK authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>
#include <memory>

void Rtrsyl3(const char *trana, const char *tranb, INTEGER const isgn, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *c, INTEGER const ldc, REAL &scale, INTEGER *iwork, INTEGER const liwork, REAL *swork, INTEGER const ldswork, INTEGER &info) {
    //
    // Decode and Test input parameters
    //
    bool notrna = Mlsame(trana, "N");
    bool notrnb = Mlsame(tranb, "N");
    //
    // Use the same block size for all matrices.
    //
    INTEGER nb = max((INTEGER)8, iMlaenv(1, "Rtrsyl", "", m, n, -1, -1));
    //
    // Compute number of blocks in A and B
    //
    INTEGER nba = max((INTEGER)1, (m + nb - 1) / nb);
    INTEGER nbb = max((INTEGER)1, (n + nb - 1) / nb);
    //
    // Compute workspace
    //
    info = 0;
    bool lquery = (liwork == -1 || ldswork == -1);
    iwork[1 - 1] = nba + nbb + 2;
    if (lquery) {
        swork[0] = max(nba, nbb);
        swork[(2 - 1)] = 2 * nbb + nba;
    }
    //
    // Test the input arguments
    //
    if (!notrna && !Mlsame(trana, "T") && !Mlsame(trana, "C")) {
        info = -1;
    } else if (!notrnb && !Mlsame(tranb, "T") && !Mlsame(tranb, "C")) {
        info = -2;
    } else if (isgn != 1 && isgn != -1) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, m)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -11;
    }
    if (info != 0) {
        Mxerbla("Rtrsyl3", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    // Quick return if possible
    //
    const REAL one = 1.0;
    const REAL two = 2.0;
    scale = one;
    if (m == 0 || n == 0) {
        return;
    }
    //
    // Use unblocked code for small problems or if insufficient
    // workspaces are provided
    //
    if (min(nba, nbb) == 1 || ldswork < max(nba, nbb) || liwork < iwork[1 - 1]) {
        Rtrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
        return;
    }
    //
    // Set constants to control overflow
    //
    REAL smlnum = Rlamch("S");
    REAL bignum = one / smlnum;
    //
    // Partition A such that 2-by-2 blocks on the diagonal are not split
    //
    bool skip = false;
    INTEGER i = 0;
    for (i = 1; i <= nba; i = i + 1) {
        iwork[i - 1] = (i - 1) * nb + 1;
    }
    iwork[(nba + 1) - 1] = m + 1;
    INTEGER k = 0;
    INTEGER l1 = 0;
    INTEGER l2 = 0;
    INTEGER l = 0;
    const REAL zero = 0.0;
    for (k = 1; k <= nba; k = k + 1) {
        l1 = iwork[k - 1];
        l2 = iwork[(k + 1) - 1] - 1;
        for (l = l1; l <= l2; l = l + 1) {
            if (skip) {
                skip = false;
                continue;
            }
            if (l >= m) {
                // A( M, M ) is a 1-by-1 block
                continue;
            }
            if (a[(l - 1) + ((l + 1) - 1) * lda] != zero && a[((l + 1) - 1) + (l - 1) * lda] != zero) {
                // Check if 2-by-2 block is split
                if (l + 1 == iwork[(k + 1) - 1]) {
                    iwork[(k + 1) - 1]++;
                    continue;
                }
                skip = true;
            }
        }
    }
    iwork[(nba + 1) - 1] = m + 1;
    if (iwork[nba - 1] >= iwork[(nba + 1) - 1]) {
        iwork[nba - 1] = iwork[(nba + 1) - 1];
        nba = nba - 1;
    }
    //
    // Partition B such that 2-by-2 blocks on the diagonal are not split
    //
    INTEGER pc = nba + 1;
    skip = false;
    for (i = 1; i <= nbb; i = i + 1) {
        iwork[(pc + i) - 1] = (i - 1) * nb + 1;
    }
    iwork[(pc + nbb + 1) - 1] = n + 1;
    for (k = 1; k <= nbb; k = k + 1) {
        l1 = iwork[(pc + k) - 1];
        l2 = iwork[(pc + k + 1) - 1] - 1;
        for (l = l1; l <= l2; l = l + 1) {
            if (skip) {
                skip = false;
                continue;
            }
            if (l >= n) {
                // B( N, N ) is a 1-by-1 block
                continue;
            }
            if (b[(l - 1) + ((l + 1) - 1) * ldb] != zero && b[((l + 1) - 1) + (l - 1) * ldb] != zero) {
                // Check if 2-by-2 block is split
                if (l + 1 == iwork[(pc + k + 1) - 1]) {
                    iwork[(pc + k + 1) - 1]++;
                    continue;
                }
                skip = true;
            }
        }
    }
    iwork[(pc + nbb + 1) - 1] = n + 1;
    if (iwork[(pc + nbb) - 1] >= iwork[(pc + nbb + 1) - 1]) {
        iwork[(pc + nbb) - 1] = iwork[(pc + nbb + 1) - 1];
        nbb = nbb - 1;
    }
    //
    // Set local scaling factors - must never attain zero.
    //
    for (l = 1; l <= nbb; l = l + 1) {
        for (k = 1; k <= nba; k = k + 1) {
            swork[(k - 1) + (l - 1) * ldswork] = one;
        }
    }
    //
    // Fallback scaling factor to prevent flushing of SWORK( K, L ) to zero.
    // This scaling is to ensure compatibility with TRSYL and may get flushed.
    //
    REAL buf = one;
    //
    // Compute upper bounds of blocks of A and B
    //
    INTEGER awrk = nbb;
    INTEGER k1 = 0;
    INTEGER k2 = 0;
    auto wnrm_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, max(m, n)));
    REAL *wnrm = wnrm_storage.get();
    for (k = 1; k <= nba; k = k + 1) {
        k1 = iwork[k - 1];
        k2 = iwork[(k + 1) - 1];
        for (l = k; l <= nba; l = l + 1) {
            l1 = iwork[l - 1];
            l2 = iwork[(l + 1) - 1];
            if (notrna) {
                swork[(k - 1) + ((awrk + l) - 1) * ldswork] = Rlange("I", k2 - k1, l2 - l1, &a[(k1 - 1) + (l1 - 1) * lda], lda, wnrm);
            } else {
                swork[(l - 1) + ((awrk + k) - 1) * ldswork] = Rlange("1", k2 - k1, l2 - l1, &a[(k1 - 1) + (l1 - 1) * lda], lda, wnrm);
            }
        }
    }
    INTEGER bwrk = nbb + nba;
    for (k = 1; k <= nbb; k = k + 1) {
        k1 = iwork[(pc + k) - 1];
        k2 = iwork[(pc + k + 1) - 1];
        for (l = k; l <= nbb; l = l + 1) {
            l1 = iwork[(pc + l) - 1];
            l2 = iwork[(pc + l + 1) - 1];
            if (notrnb) {
                swork[(k - 1) + ((bwrk + l) - 1) * ldswork] = Rlange("I", k2 - k1, l2 - l1, &b[(k1 - 1) + (l1 - 1) * ldb], ldb, wnrm);
            } else {
                swork[(l - 1) + ((bwrk + k) - 1) * ldswork] = Rlange("1", k2 - k1, l2 - l1, &b[(k1 - 1) + (l1 - 1) * ldb], ldb, wnrm);
            }
        }
    }
    //
    REAL sgn = castREAL(isgn);
    //
    REAL scaloc = 0.0;
    INTEGER iinfo = 0;
    INTEGER jj = 0;
    INTEGER ll = 0;
    REAL xnrm = 0.0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    REAL cnrm = 0.0;
    REAL scamin = 0.0;
    REAL anrm = 0.0;
    REAL scal = 0.0;
    INTEGER j = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    REAL bnrm = 0.0;
    if (notrna && notrnb) {
        //
        // Solve    A*X + ISGN*X*B = scale*C.
        //
        // The (K,L)th block of X is determined starting from
        // bottom-left corner column by column by
        //
        // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
        //
        // Where
        // M                         L-1
        // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
        // I=K+1                       J=1
        //
        // Start loop over block rows (index = K) and block columns (index = L)
        //
        for (k = nba; k >= 1; k = k - 1) {
            //
            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )
            //
            k1 = iwork[k - 1];
            k2 = iwork[(k + 1) - 1];
            for (l = 1; l <= nbb; l = l + 1) {
                //
                // L1: column index of the first column in X( K, L )
                // L2: column index of the first column in X( K, L + 1)
                // so that L2 - L1 is the row count of the block X( K, L )
                //
                l1 = iwork[(pc + l) - 1];
                l2 = iwork[(pc + l + 1) - 1];
                //
                Rtrsyl(trana, tranb, isgn, k2 - k1, l2 - l1, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, scaloc, iinfo);
                info = max(info, iinfo);
                //
                if (scaloc * swork[(k - 1) + (l - 1) * ldswork] == zero) {
                    if (scaloc == zero) {
                        // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                        // is larger than the product of BIGNUM**2 and cannot be
                        // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                        // Mark the computation as pointless.
                        buf = zero;
                    } else {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                    }
                    for (jj = 1; jj <= nbb; jj = jj + 1) {
                        for (ll = 1; ll <= nba; ll = ll + 1) {
                            // Bound by BIGNUM to not introduce Inf. The value
                            // is irrelevant; corresponding entries of the
                            // solution will be flushed in consistency scaling.
                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                        }
                    }
                }
                swork[(k - 1) + (l - 1) * ldswork] = scaloc * swork[(k - 1) + (l - 1) * ldswork];
                xnrm = Rlange("I", k2 - k1, l2 - l1, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                //
                for (i = k - 1; i >= 1; i = i - 1) {
                    //
                    // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )
                    //
                    i1 = iwork[i - 1];
                    i2 = iwork[(i + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", i2 - i1, l2 - l1, &c[(i1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(i - 1) + (l - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(i - 1) + (l - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    anrm = swork[(i - 1) + ((awrk + k) - 1) * ldswork];
                    scaloc = Rlarmm(anrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to C( I, L ) and C( K, L ).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (jj = l1; jj <= l2 - 1; jj = jj + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (jj - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(i - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(i2 - i1, scal, &c[(i1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(i - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("N", "N", i2 - i1, l2 - l1, k2 - k1, -one, &a[(i1 - 1) + (k1 - 1) * lda], lda, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, one, &c[(i1 - 1) + (l1 - 1) * ldc], ldc);
                    //
                }
                //
                for (j = l + 1; j <= nbb; j = j + 1) {
                    //
                    // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )
                    //
                    j1 = iwork[(pc + j) - 1];
                    j2 = iwork[(pc + j + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", k2 - k1, j2 - j1, &c[(k1 - 1) + (j1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(k - 1) + (j - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(k - 1) + (j - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    bnrm = swork[(l - 1) + ((bwrk + j) - 1) * ldswork];
                    scaloc = Rlarmm(bnrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to C( K, J ) and C( K, L).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(k - 1) + (j - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (jj = j1; jj <= j2 - 1; jj = jj + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (jj - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(k - 1) + (j - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("N", "N", k2 - k1, j2 - j1, l2 - l1, -sgn, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, &b[(l1 - 1) + (j1 - 1) * ldb], ldb, one, &c[(k1 - 1) + (j1 - 1) * ldc], ldc);
                }
            }
        }
    } else if (!notrna && notrnb) {
        //
        // Solve    A**T*X + ISGN*X*B = scale*C.
        //
        // The (K,L)th block of X is determined starting from
        // upper-left corner column by column by
        //
        // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
        //
        // Where
        // K-1                        L-1
        // R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
        // I=1                        J=1
        //
        // Start loop over block rows (index = K) and block columns (index = L)
        //
        for (k = 1; k <= nba; k = k + 1) {
            //
            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )
            //
            k1 = iwork[k - 1];
            k2 = iwork[(k + 1) - 1];
            for (l = 1; l <= nbb; l = l + 1) {
                //
                // L1: column index of the first column in X( K, L )
                // L2: column index of the first column in X( K, L + 1)
                // so that L2 - L1 is the row count of the block X( K, L )
                //
                l1 = iwork[(pc + l) - 1];
                l2 = iwork[(pc + l + 1) - 1];
                //
                Rtrsyl(trana, tranb, isgn, k2 - k1, l2 - l1, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, scaloc, iinfo);
                info = max(info, iinfo);
                //
                if (scaloc * swork[(k - 1) + (l - 1) * ldswork] == zero) {
                    if (scaloc == zero) {
                        // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                        // is larger than the product of BIGNUM**2 and cannot be
                        // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                        // Mark the computation as pointless.
                        buf = zero;
                    } else {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                    }
                    for (jj = 1; jj <= nbb; jj = jj + 1) {
                        for (ll = 1; ll <= nba; ll = ll + 1) {
                            // Bound by BIGNUM to not introduce Inf. The value
                            // is irrelevant; corresponding entries of the
                            // solution will be flushed in consistency scaling.
                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                        }
                    }
                }
                swork[(k - 1) + (l - 1) * ldswork] = scaloc * swork[(k - 1) + (l - 1) * ldswork];
                xnrm = Rlange("I", k2 - k1, l2 - l1, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                //
                for (i = k + 1; i <= nba; i = i + 1) {
                    //
                    // C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L )
                    //
                    i1 = iwork[i - 1];
                    i2 = iwork[(i + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", i2 - i1, l2 - l1, &c[(i1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(i - 1) + (l - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(i - 1) + (l - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    anrm = swork[(i - 1) + ((awrk + k) - 1) * ldswork];
                    scaloc = Rlarmm(anrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to to C( I, L ) and C( K, L ).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(i - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(i2 - i1, scal, &c[(i1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(i - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("T", "N", i2 - i1, l2 - l1, k2 - k1, -one, &a[(k1 - 1) + (i1 - 1) * lda], lda, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, one, &c[(i1 - 1) + (l1 - 1) * ldc], ldc);
                }
                //
                for (j = l + 1; j <= nbb; j = j + 1) {
                    //
                    // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )
                    //
                    j1 = iwork[(pc + j) - 1];
                    j2 = iwork[(pc + j + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", k2 - k1, j2 - j1, &c[(k1 - 1) + (j1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(k - 1) + (j - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(k - 1) + (j - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    bnrm = swork[(l - 1) + ((bwrk + j) - 1) * ldswork];
                    scaloc = Rlarmm(bnrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to to C( K, J ) and C( K, L ).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(k - 1) + (j - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (jj = j1; jj <= j2 - 1; jj = jj + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (jj - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(k - 1) + (j - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("N", "N", k2 - k1, j2 - j1, l2 - l1, -sgn, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, &b[(l1 - 1) + (j1 - 1) * ldb], ldb, one, &c[(k1 - 1) + (j1 - 1) * ldc], ldc);
                }
            }
        }
    } else if (!notrna && !notrnb) {
        //
        // Solve    A**T*X + ISGN*X*B**T = scale*C.
        //
        // The (K,L)th block of X is determined starting from
        // top-right corner column by column by
        //
        // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
        //
        // Where
        // K-1                          N
        // R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
        // I=1                        J=L+1
        //
        // Start loop over block rows (index = K) and block columns (index = L)
        //
        for (k = 1; k <= nba; k = k + 1) {
            //
            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )
            //
            k1 = iwork[k - 1];
            k2 = iwork[(k + 1) - 1];
            for (l = nbb; l >= 1; l = l - 1) {
                //
                // L1: column index of the first column in X( K, L )
                // L2: column index of the first column in X( K, L + 1)
                // so that L2 - L1 is the row count of the block X( K, L )
                //
                l1 = iwork[(pc + l) - 1];
                l2 = iwork[(pc + l + 1) - 1];
                //
                Rtrsyl(trana, tranb, isgn, k2 - k1, l2 - l1, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, scaloc, iinfo);
                info = max(info, iinfo);
                //
                swork[(k - 1) + (l - 1) * ldswork] = scaloc * swork[(k - 1) + (l - 1) * ldswork];
                if (scaloc * swork[(k - 1) + (l - 1) * ldswork] == zero) {
                    if (scaloc == zero) {
                        // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                        // is larger than the product of BIGNUM**2 and cannot be
                        // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                        // Mark the computation as pointless.
                        buf = zero;
                    } else {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                    }
                    for (jj = 1; jj <= nbb; jj = jj + 1) {
                        for (ll = 1; ll <= nba; ll = ll + 1) {
                            // Bound by BIGNUM to not introduce Inf. The value
                            // is irrelevant; corresponding entries of the
                            // solution will be flushed in consistency scaling.
                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                        }
                    }
                }
                xnrm = Rlange("I", k2 - k1, l2 - l1, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                //
                for (i = k + 1; i <= nba; i = i + 1) {
                    //
                    // C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L )
                    //
                    i1 = iwork[i - 1];
                    i2 = iwork[(i + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", i2 - i1, l2 - l1, &c[(i1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(i - 1) + (l - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(i - 1) + (l - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    anrm = swork[(i - 1) + ((awrk + k) - 1) * ldswork];
                    scaloc = Rlarmm(anrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to C( I, L ) and C( K, L ).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(i - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(i2 - i1, scal, &c[(i1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(i - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("T", "N", i2 - i1, l2 - l1, k2 - k1, -one, &a[(k1 - 1) + (i1 - 1) * lda], lda, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, one, &c[(i1 - 1) + (l1 - 1) * ldc], ldc);
                }
                //
                for (j = 1; j <= l - 1; j = j + 1) {
                    //
                    // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T
                    //
                    j1 = iwork[(pc + j) - 1];
                    j2 = iwork[(pc + j + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", k2 - k1, j2 - j1, &c[(k1 - 1) + (j1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(k - 1) + (j - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(k - 1) + (j - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    bnrm = swork[(l - 1) + ((bwrk + j) - 1) * ldswork];
                    scaloc = Rlarmm(bnrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to C( K, J ) and C( K, L ).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(k - 1) + (j - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (jj = j1; jj <= j2 - 1; jj = jj + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (jj - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(k - 1) + (j - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("N", "T", k2 - k1, j2 - j1, l2 - l1, -sgn, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, &b[(j1 - 1) + (l1 - 1) * ldb], ldb, one, &c[(k1 - 1) + (j1 - 1) * ldc], ldc);
                }
            }
        }
    } else if (notrna && !notrnb) {
        //
        // Solve    A*X + ISGN*X*B**T = scale*C.
        //
        // The (K,L)th block of X is determined starting from
        // bottom-right corner column by column by
        //
        // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
        //
        // Where
        // M                          N
        // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
        // I=K+1                      J=L+1
        //
        // Start loop over block rows (index = K) and block columns (index = L)
        //
        for (k = nba; k >= 1; k = k - 1) {
            //
            // K1: row index of the first row in X( K, L )
            // K2: row index of the first row in X( K+1, L )
            // so the K2 - K1 is the column count of the block X( K, L )
            //
            k1 = iwork[k - 1];
            k2 = iwork[(k + 1) - 1];
            for (l = nbb; l >= 1; l = l - 1) {
                //
                // L1: column index of the first column in X( K, L )
                // L2: column index of the first column in X( K, L + 1)
                // so that L2 - L1 is the row count of the block X( K, L )
                //
                l1 = iwork[(pc + l) - 1];
                l2 = iwork[(pc + l + 1) - 1];
                //
                Rtrsyl(trana, tranb, isgn, k2 - k1, l2 - l1, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, scaloc, iinfo);
                info = max(info, iinfo);
                //
                if (scaloc * swork[(k - 1) + (l - 1) * ldswork] == zero) {
                    if (scaloc == zero) {
                        // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
                        // is larger than the product of BIGNUM**2 and cannot be
                        // represented in the form (1/SCALE)*X(K1:K2-1, L1:L2-1).
                        // Mark the computation as pointless.
                        buf = zero;
                    } else {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                    }
                    for (jj = 1; jj <= nbb; jj = jj + 1) {
                        for (ll = 1; ll <= nba; ll = ll + 1) {
                            // Bound by BIGNUM to not introduce Inf. The value
                            // is irrelevant; corresponding entries of the
                            // solution will be flushed in consistency scaling.
                            swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                        }
                    }
                }
                swork[(k - 1) + (l - 1) * ldswork] = scaloc * swork[(k - 1) + (l - 1) * ldswork];
                xnrm = Rlange("I", k2 - k1, l2 - l1, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                //
                for (i = 1; i <= k - 1; i = i + 1) {
                    //
                    // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )
                    //
                    i1 = iwork[i - 1];
                    i2 = iwork[(i + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", i2 - i1, l2 - l1, &c[(i1 - 1) + (l1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(i - 1) + (l - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(i - 1) + (l - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    anrm = swork[(i - 1) + ((awrk + k) - 1) * ldswork];
                    scaloc = Rlarmm(anrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to C( I, L ) and C( K, L ).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(i - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                            Rscal(i2 - i1, scal, &c[(i1 - 1) + (ll - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(i - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("N", "N", i2 - i1, l2 - l1, k2 - k1, -one, &a[(i1 - 1) + (k1 - 1) * lda], lda, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, one, &c[(i1 - 1) + (l1 - 1) * ldc], ldc);
                    //
                }
                //
                for (j = 1; j <= l - 1; j = j + 1) {
                    //
                    // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T
                    //
                    j1 = iwork[(pc + j) - 1];
                    j2 = iwork[(pc + j + 1) - 1];
                    //
                    // Compute scaling factor to survive the linear update
                    // simulating consistent scaling.
                    //
                    cnrm = Rlange("I", k2 - k1, j2 - j1, &c[(k1 - 1) + (j1 - 1) * ldc], ldc, wnrm);
                    scamin = min(swork[(k - 1) + (j - 1) * ldswork], swork[(k - 1) + (l - 1) * ldswork]);
                    cnrm = cnrm * (scamin / swork[(k - 1) + (j - 1) * ldswork]);
                    xnrm = xnrm * (scamin / swork[(k - 1) + (l - 1) * ldswork]);
                    bnrm = swork[(l - 1) + ((bwrk + j) - 1) * ldswork];
                    scaloc = Rlarmm(bnrm, xnrm, cnrm);
                    if (scaloc * scamin == zero) {
                        // Use second scaling factor to prevent flushing to zero.
                        buf = buf * pow(two, Mexponent(scaloc));
                        for (jj = 1; jj <= nbb; jj = jj + 1) {
                            for (ll = 1; ll <= nba; ll = ll + 1) {
                                swork[(ll - 1) + (jj - 1) * ldswork] = min(bignum, swork[(ll - 1) + (jj - 1) * ldswork] / pow(two, Mexponent(scaloc)));
                            }
                        }
                        scamin = scamin / pow(two, Mexponent(scaloc));
                        scaloc = scaloc / pow(two, Mexponent(scaloc));
                    }
                    cnrm = cnrm * scaloc;
                    xnrm = xnrm * scaloc;
                    //
                    // Simultaneously apply the robust update factor and the
                    // consistency scaling factor to C( K, J ) and C( K, L ).
                    //
                    scal = (scamin / swork[(k - 1) + (l - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (jj = l1; jj <= l2 - 1; jj = jj + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (jj - 1) * ldc], 1);
                        }
                    }
                    //
                    scal = (scamin / swork[(k - 1) + (j - 1) * ldswork]) * scaloc;
                    if (scal != one) {
                        for (jj = j1; jj <= j2 - 1; jj = jj + 1) {
                            Rscal(k2 - k1, scal, &c[(k1 - 1) + (jj - 1) * ldc], 1);
                        }
                    }
                    //
                    // Record current scaling factor
                    //
                    swork[(k - 1) + (l - 1) * ldswork] = scamin * scaloc;
                    swork[(k - 1) + (j - 1) * ldswork] = scamin * scaloc;
                    //
                    Rgemm("N", "T", k2 - k1, j2 - j1, l2 - l1, -sgn, &c[(k1 - 1) + (l1 - 1) * ldc], ldc, &b[(j1 - 1) + (l1 - 1) * ldb], ldb, one, &c[(k1 - 1) + (j1 - 1) * ldc], ldc);
                }
            }
        }
        //
    }
    //
    // Reduce local scaling factors
    //
    scale = swork[0];
    for (k = 1; k <= nba; k = k + 1) {
        for (l = 1; l <= nbb; l = l + 1) {
            scale = min(scale, swork[(k - 1) + (l - 1) * ldswork]);
        }
    }
    //
    if (scale == zero) {
        //
        // The magnitude of the largest entry of the solution is larger
        // than the product of BIGNUM**2 and cannot be represented in the
        // form (1/SCALE)*X if SCALE is DOUBLE PRECISION. Set SCALE to
        // zero and give up.
        //
        iwork[1 - 1] = nba + nbb + 2;
        swork[0] = max(nba, nbb);
        swork[(2 - 1)] = 2 * nbb + nba;
        return;
    }
    //
    // Realize consistent scaling
    //
    for (k = 1; k <= nba; k = k + 1) {
        k1 = iwork[k - 1];
        k2 = iwork[(k + 1) - 1];
        for (l = 1; l <= nbb; l = l + 1) {
            l1 = iwork[(pc + l) - 1];
            l2 = iwork[(pc + l + 1) - 1];
            scal = scale / swork[(k - 1) + (l - 1) * ldswork];
            if (scal != one) {
                for (ll = l1; ll <= l2 - 1; ll = ll + 1) {
                    Rscal(k2 - k1, scal, &c[(k1 - 1) + (ll - 1) * ldc], 1);
                }
            }
        }
    }
    //
    if (buf != one && buf > zero) {
        //
        // Decrease SCALE as much as possible.
        //
        scaloc = min(scale / smlnum, one / buf);
        buf = buf * scaloc;
        scale = scale / scaloc;
    }
    //
    if (buf != one && buf > zero) {
        //
        // In case of overly aggressive scaling during the computation,
        // flushing of the global scale factor may be prevented by
        // undoing some of the scaling. This step is to ensure that
        // this routine flushes only scale factors that TRSYL also
        // flushes and be usable as a drop-in replacement.
        //
        // How much can the normwise largest entry be upscaled?
        //
        scal = c[0];
        for (k = 1; k <= m; k = k + 1) {
            for (l = 1; l <= n; l = l + 1) {
                scal = max(scal, abs(c[(k - 1) + (l - 1) * ldc]));
            }
        }
        //
        // Increase BUF as close to 1 as possible and apply scaling.
        //
        scaloc = min(bignum / scal, one / buf);
        buf = buf * scaloc;
        Rlascl("G", -1, -1, one, scaloc, m, n, c, ldc, iwork[1 - 1]);
    }
    //
    // Combine with buffer scaling factor. SCALE will be flushed if
    // BUF is less than one here.
    //
    scale = scale * buf;
    //
    // Restore workspace dimensions
    //
    iwork[1 - 1] = nba + nbb + 2;
    swork[0] = max(nba, nbb);
    swork[(2 - 1)] = 2 * nbb + nba;
    //
    // End of Rtrsyl3
    //
}

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
#include <mplapack_matgen.h>

void Clatmr(INTEGER const m, INTEGER const n, const char *dist, INTEGER *iseed, const char *sym, COMPLEX *d, INTEGER const mode, REAL const cond, COMPLEX const dmax, const char *rsign, const char *grade, COMPLEX *dl, INTEGER const model, REAL const condl, COMPLEX *dr, INTEGER const moder, REAL const condr, const char *pivtng, INTEGER *ipivot, INTEGER const kl, INTEGER const ku, REAL const sparse, REAL const anorm, const char *pack, COMPLEX *a, INTEGER const lda, INTEGER *iwork, INTEGER &info) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     1)      Decode and Test the input parameters.
    //             Initialize flags & seed.
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Decode DIST
    //
    INTEGER idist = 0;
    if (Mlsame(dist, "U")) {
        idist = 1;
    } else if (Mlsame(dist, "S")) {
        idist = 2;
    } else if (Mlsame(dist, "N")) {
        idist = 3;
    } else if (Mlsame(dist, "D")) {
        idist = 4;
    } else {
        idist = -1;
    }
    //
    //     Decode SYM
    //
    INTEGER isym = 0;
    if (Mlsame(sym, "H")) {
        isym = 0;
    } else if (Mlsame(sym, "N")) {
        isym = 1;
    } else if (Mlsame(sym, "S")) {
        isym = 2;
    } else {
        isym = -1;
    }
    //
    //     Decode RSIGN
    //
    INTEGER irsign = 0;
    if (Mlsame(rsign, "F")) {
        irsign = 0;
    } else if (Mlsame(rsign, "T")) {
        irsign = 1;
    } else {
        irsign = -1;
    }
    //
    //     Decode PIVTNG
    //
    INTEGER ipvtng = 0;
    INTEGER npvts = 0;
    if (Mlsame(pivtng, "N")) {
        ipvtng = 0;
    } else if (Mlsame(pivtng, " ")) {
        ipvtng = 0;
    } else if (Mlsame(pivtng, "L")) {
        ipvtng = 1;
        npvts = m;
    } else if (Mlsame(pivtng, "R")) {
        ipvtng = 2;
        npvts = n;
    } else if (Mlsame(pivtng, "B")) {
        ipvtng = 3;
        npvts = min(n, m);
    } else if (Mlsame(pivtng, "F")) {
        ipvtng = 3;
        npvts = min(n, m);
    } else {
        ipvtng = -1;
    }
    //
    //     Decode GRADE
    //
    INTEGER igrade = 0;
    if (Mlsame(grade, "N")) {
        igrade = 0;
    } else if (Mlsame(grade, "L")) {
        igrade = 1;
    } else if (Mlsame(grade, "R")) {
        igrade = 2;
    } else if (Mlsame(grade, "B")) {
        igrade = 3;
    } else if (Mlsame(grade, "E")) {
        igrade = 4;
    } else if (Mlsame(grade, "H")) {
        igrade = 5;
    } else if (Mlsame(grade, "S")) {
        igrade = 6;
    } else {
        igrade = -1;
    }
    //
    //     Decode PACK
    //
    INTEGER ipack = 0;
    if (Mlsame(pack, "N")) {
        ipack = 0;
    } else if (Mlsame(pack, "U")) {
        ipack = 1;
    } else if (Mlsame(pack, "L")) {
        ipack = 2;
    } else if (Mlsame(pack, "C")) {
        ipack = 3;
    } else if (Mlsame(pack, "R")) {
        ipack = 4;
    } else if (Mlsame(pack, "B")) {
        ipack = 5;
    } else if (Mlsame(pack, "Q")) {
        ipack = 6;
    } else if (Mlsame(pack, "Z")) {
        ipack = 7;
    } else {
        ipack = -1;
    }
    //
    //     Set certain internal parameters
    //
    INTEGER mnmin = min(m, n);
    INTEGER kll = min(kl, m - 1);
    INTEGER kuu = min(ku, n - 1);
    //
    //     If inv(DL) is used, check to see if DL has a zero entry.
    //
    bool dzero = false;
    INTEGER i = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (igrade == 4 && model == 0) {
        for (i = 1; i <= m; i = i + 1) {
            if (dl[i - 1] == czero) {
                dzero = true;
            }
        }
    }
    //
    //     Check values in IPIVOT
    //
    bool badpvt = false;
    INTEGER j = 0;
    if (ipvtng > 0) {
        for (j = 1; j <= npvts; j = j + 1) {
            if (ipivot[j - 1] <= 0 || ipivot[j - 1] > npvts) {
                badpvt = true;
            }
        }
    }
    //
    //     Set INFO if an error
    //
    const REAL one = 1.0;
    const REAL zero = 0.0;
    if (m < 0) {
        info = -1;
    } else if (m != n && (isym == 0 || isym == 2)) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (idist == -1) {
        info = -3;
    } else if (isym == -1) {
        info = -5;
    } else if (mode < -6 || mode > 6) {
        info = -7;
    } else if ((mode != -6 && mode != 0 && mode != 6) && cond < one) {
        info = -8;
    } else if ((mode != -6 && mode != 0 && mode != 6) && irsign == -1) {
        info = -10;
    } else if (igrade == -1 || (igrade == 4 && m != n) || ((igrade == 1 || igrade == 2 || igrade == 3 || igrade == 4 || igrade == 6) && isym == 0) || ((igrade == 1 || igrade == 2 || igrade == 3 || igrade == 4 || igrade == 5) && isym == 2)) {
        info = -11;
    } else if (igrade == 4 && dzero) {
        info = -12;
    } else if ((igrade == 1 || igrade == 3 || igrade == 4 || igrade == 5 || igrade == 6) && (model < -6 || model > 6)) {
        info = -13;
    } else if ((igrade == 1 || igrade == 3 || igrade == 4 || igrade == 5 || igrade == 6) && (model != -6 && model != 0 && model != 6) && condl < one) {
        info = -14;
    } else if ((igrade == 2 || igrade == 3) && (moder < -6 || moder > 6)) {
        info = -16;
    } else if ((igrade == 2 || igrade == 3) && (moder != -6 && moder != 0 && moder != 6) && condr < one) {
        info = -17;
    } else if (ipvtng == -1 || (ipvtng == 3 && m != n) || ((ipvtng == 1 || ipvtng == 2) && (isym == 0 || isym == 2))) {
        info = -18;
    } else if (ipvtng != 0 && badpvt) {
        info = -19;
    } else if (kl < 0) {
        info = -20;
    } else if (ku < 0 || ((isym == 0 || isym == 2) && kl != ku)) {
        info = -21;
    } else if (sparse < zero || sparse > one) {
        info = -22;
    } else if (ipack == -1 || ((ipack == 1 || ipack == 2 || ipack == 5 || ipack == 6) && isym == 1) || (ipack == 3 && isym == 1 && (kl != 0 || m != n)) || (ipack == 4 && isym == 1 && (ku != 0 || m != n))) {
        info = -24;
    } else if (((ipack == 0 || ipack == 1 || ipack == 2) && lda < max((INTEGER)1, m)) || ((ipack == 3 || ipack == 4) && lda < 1) || ((ipack == 5 || ipack == 6) && lda < kuu + 1) || (ipack == 7 && lda < kll + kuu + 1)) {
        info = -26;
    }
    //
    if (info != 0) {
        Mxerbla("Clatmr", -info);
        return;
    }
    //
    //     Decide if we can pivot consistently
    //
    bool fulbnd = false;
    if (kuu == n - 1 && kll == m - 1) {
        fulbnd = true;
    }
    //
    //     Initialize random number generator
    //
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = mod(abs(iseed[i - 1]), 4096);
    }
    //
    iseed[4 - 1] = 2 * (iseed[4 - 1] / 2) + 1;
    //
    //     2)      Set up D, DL, and DR, if indicated.
    //
    //             Compute D according to COND and MODE
    //
    Clatm1(mode, cond, irsign, idist, iseed, d, mnmin, info);
    if (info != 0) {
        info = 1;
        return;
    }
    REAL temp = 0.0;
    COMPLEX calpha = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    if (mode != 0 && mode != -6 && mode != 6) {
        //
        //        Scale by DMAX
        //
        temp = abs(d[1 - 1]);
        for (i = 2; i <= mnmin; i = i + 1) {
            temp = max(temp, abs(d[i - 1]));
        }
        if (temp == zero && dmax != czero) {
            info = 2;
            return;
        }
        if (temp != zero) {
            calpha = dmax / temp;
        } else {
            calpha = cone;
        }
        for (i = 1; i <= mnmin; i = i + 1) {
            d[i - 1] = calpha * d[i - 1];
        }
        //
    }
    //
    //     If matrix Hermitian, make D real
    //
    if (isym == 0) {
        for (i = 1; i <= mnmin; i = i + 1) {
            d[i - 1] = d[i - 1].real();
        }
    }
    //
    //     Compute DL if grading set
    //
    if (igrade == 1 || igrade == 3 || igrade == 4 || igrade == 5 || igrade == 6) {
        Clatm1(model, condl, 0, idist, iseed, dl, m, info);
        if (info != 0) {
            info = 3;
            return;
        }
    }
    //
    //     Compute DR if grading set
    //
    if (igrade == 2 || igrade == 3) {
        Clatm1(moder, condr, 0, idist, iseed, dr, n, info);
        if (info != 0) {
            info = 4;
            return;
        }
    }
    //
    //     3)     Generate IWORK if pivoting
    //
    INTEGER k = 0;
    if (ipvtng > 0) {
        for (i = 1; i <= npvts; i = i + 1) {
            iwork[i - 1] = i;
        }
        if (fulbnd) {
            for (i = 1; i <= npvts; i = i + 1) {
                k = ipivot[i - 1];
                j = iwork[i - 1];
                iwork[i - 1] = iwork[k - 1];
                iwork[k - 1] = j;
            }
        } else {
            for (i = npvts; i >= 1; i = i - 1) {
                k = ipivot[i - 1];
                j = iwork[i - 1];
                iwork[i - 1] = iwork[k - 1];
                iwork[k - 1] = j;
            }
        }
    }
    //
    //     4)      Generate matrices for each kind of PACKing
    //             Always sweep matrix columnwise (if symmetric, upper
    //             half only) so that matrix generated does not depend
    //             on PACK
    //
    INTEGER isub = 0;
    INTEGER jsub = 0;
    COMPLEX ctemp = 0.0;
    INTEGER mnsub = 0;
    INTEGER mxsub = 0;
    INTEGER jjsub = 0;
    INTEGER iisub = 0;
    if (fulbnd) {
        //
        //        Use Clatm3 so matrices generated with differing PIVOTing only
        //        differ only in the order of their rows and/or columns.
        //
        if (ipack == 0) {
            if (isym == 0) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        a[(isub - 1) + (jsub - 1) * lda] = ctemp;
                        a[(jsub - 1) + (isub - 1) * lda] = conj(ctemp);
                    }
                }
            } else if (isym == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= m; i = i + 1) {
                        ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        a[(isub - 1) + (jsub - 1) * lda] = ctemp;
                    }
                }
            } else if (isym == 2) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        a[(isub - 1) + (jsub - 1) * lda] = ctemp;
                        a[(jsub - 1) + (isub - 1) * lda] = ctemp;
                    }
                }
            }
            //
        } else if (ipack == 1) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    mnsub = min(isub, jsub);
                    mxsub = max(isub, jsub);
                    if (mxsub == isub && isym == 0) {
                        a[(mnsub - 1) + (mxsub - 1) * lda] = conj(ctemp);
                    } else {
                        a[(mnsub - 1) + (mxsub - 1) * lda] = ctemp;
                    }
                    if (mnsub != mxsub) {
                        a[(mxsub - 1) + (mnsub - 1) * lda] = czero;
                    }
                }
            }
            //
        } else if (ipack == 2) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    mnsub = min(isub, jsub);
                    mxsub = max(isub, jsub);
                    if (mxsub == jsub && isym == 0) {
                        a[(mxsub - 1) + (mnsub - 1) * lda] = conj(ctemp);
                    } else {
                        a[(mxsub - 1) + (mnsub - 1) * lda] = ctemp;
                    }
                    if (mnsub != mxsub) {
                        a[(mnsub - 1) + (mxsub - 1) * lda] = czero;
                    }
                }
            }
            //
        } else if (ipack == 3) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    //
                    //                 Compute K = location of (ISUB,JSUB) entry in packed
                    //                 array
                    //
                    mnsub = min(isub, jsub);
                    mxsub = max(isub, jsub);
                    k = mxsub * (mxsub - 1) / 2 + mnsub;
                    //
                    //                 Convert K to (IISUB,JJSUB) location
                    //
                    jjsub = (k - 1) / lda + 1;
                    iisub = k - lda * (jjsub - 1);
                    //
                    if (mxsub == isub && isym == 0) {
                        a[(iisub - 1) + (jjsub - 1) * lda] = conj(ctemp);
                    } else {
                        a[(iisub - 1) + (jjsub - 1) * lda] = ctemp;
                    }
                }
            }
            //
        } else if (ipack == 4) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    //
                    //                 Compute K = location of (I,J) entry in packed array
                    //
                    mnsub = min(isub, jsub);
                    mxsub = max(isub, jsub);
                    if (mnsub == 1) {
                        k = mxsub;
                    } else {
                        k = n * (n + 1) / 2 - (n - mnsub + 1) * (n - mnsub + 2) / 2 + mxsub - mnsub + 1;
                    }
                    //
                    //                 Convert K to (IISUB,JJSUB) location
                    //
                    jjsub = (k - 1) / lda + 1;
                    iisub = k - lda * (jjsub - 1);
                    //
                    if (mxsub == jsub && isym == 0) {
                        a[(iisub - 1) + (jjsub - 1) * lda] = conj(ctemp);
                    } else {
                        a[(iisub - 1) + (jjsub - 1) * lda] = ctemp;
                    }
                }
            }
            //
        } else if (ipack == 5) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = j - kuu; i <= j; i = i + 1) {
                    if (i < 1) {
                        a[((j - i + 1) - 1) + ((i + n) - 1) * lda] = czero;
                    } else {
                        ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        mnsub = min(isub, jsub);
                        mxsub = max(isub, jsub);
                        if (mxsub == jsub && isym == 0) {
                            a[((mxsub - mnsub + 1) - 1) + (mnsub - 1) * lda] = conj(ctemp);
                        } else {
                            a[((mxsub - mnsub + 1) - 1) + (mnsub - 1) * lda] = ctemp;
                        }
                    }
                }
            }
            //
        } else if (ipack == 6) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = j - kuu; i <= j; i = i + 1) {
                    ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    mnsub = min(isub, jsub);
                    mxsub = max(isub, jsub);
                    if (mxsub == isub && isym == 0) {
                        a[((mnsub - mxsub + kuu + 1) - 1) + (mxsub - 1) * lda] = conj(ctemp);
                    } else {
                        a[((mnsub - mxsub + kuu + 1) - 1) + (mxsub - 1) * lda] = ctemp;
                    }
                }
            }
            //
        } else if (ipack == 7) {
            //
            if (isym != 1) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j - kuu; i <= j; i = i + 1) {
                        ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        mnsub = min(isub, jsub);
                        mxsub = max(isub, jsub);
                        if (i < 1) {
                            a[((j - i + 1 + kuu) - 1) + ((i + n) - 1) * lda] = czero;
                        }
                        if (mxsub == isub && isym == 0) {
                            a[((mnsub - mxsub + kuu + 1) - 1) + (mxsub - 1) * lda] = conj(ctemp);
                        } else {
                            a[((mnsub - mxsub + kuu + 1) - 1) + (mxsub - 1) * lda] = ctemp;
                        }
                        if (i >= 1 && mnsub != mxsub) {
                            if (mnsub == isub && isym == 0) {
                                a[((mxsub - mnsub + 1 + kuu) - 1) + (mnsub - 1) * lda] = conj(ctemp);
                            } else {
                                a[((mxsub - mnsub + 1 + kuu) - 1) + (mnsub - 1) * lda] = ctemp;
                            }
                        }
                    }
                }
            } else if (isym == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j - kuu; i <= j + kll; i = i + 1) {
                        ctemp = Clatm3(m, n, i, j, isub, jsub, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        a[((isub - jsub + kuu + 1) - 1) + (jsub - 1) * lda] = ctemp;
                    }
                }
            }
            //
        }
        //
    } else {
        //
        //        Use Clatm2
        //
        if (ipack == 0) {
            if (isym == 0) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        a[(j - 1) + (i - 1) * lda] = conj(a[(i - 1) + (j - 1) * lda]);
                    }
                }
            } else if (isym == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= m; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    }
                }
            } else if (isym == 2) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        a[(j - 1) + (i - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                    }
                }
            }
            //
        } else if (ipack == 1) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    a[(i - 1) + (j - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    if (i != j) {
                        a[(j - 1) + (i - 1) * lda] = czero;
                    }
                }
            }
            //
        } else if (ipack == 2) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    if (isym == 0) {
                        a[(j - 1) + (i - 1) * lda] = conj(Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse));
                    } else {
                        a[(j - 1) + (i - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    }
                    if (i != j) {
                        a[(i - 1) + (j - 1) * lda] = czero;
                    }
                }
            }
            //
        } else if (ipack == 3) {
            //
            isub = 0;
            jsub = 1;
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    isub++;
                    if (isub > lda) {
                        isub = 1;
                        jsub++;
                    }
                    a[(isub - 1) + (jsub - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                }
            }
            //
        } else if (ipack == 4) {
            //
            if (isym == 0 || isym == 2) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        //
                        //                    Compute K = location of (I,J) entry in packed array
                        //
                        if (i == 1) {
                            k = j;
                        } else {
                            k = n * (n + 1) / 2 - (n - i + 1) * (n - i + 2) / 2 + j - i + 1;
                        }
                        //
                        //                    Convert K to (ISUB,JSUB) location
                        //
                        jsub = (k - 1) / lda + 1;
                        isub = k - lda * (jsub - 1);
                        //
                        a[(isub - 1) + (jsub - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        if (isym == 0) {
                            a[(isub - 1) + (jsub - 1) * lda] = conj(a[(isub - 1) + (jsub - 1) * lda]);
                        }
                    }
                }
            } else {
                isub = 0;
                jsub = 1;
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j; i <= m; i = i + 1) {
                        isub++;
                        if (isub > lda) {
                            isub = 1;
                            jsub++;
                        }
                        a[(isub - 1) + (jsub - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    }
                }
            }
            //
        } else if (ipack == 5) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = j - kuu; i <= j; i = i + 1) {
                    if (i < 1) {
                        a[((j - i + 1) - 1) + ((i + n) - 1) * lda] = czero;
                    } else {
                        if (isym == 0) {
                            a[((j - i + 1) - 1) + (i - 1) * lda] = conj(Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse));
                        } else {
                            a[((j - i + 1) - 1) + (i - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        }
                    }
                }
            }
            //
        } else if (ipack == 6) {
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = j - kuu; i <= j; i = i + 1) {
                    a[((i - j + kuu + 1) - 1) + (j - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                }
            }
            //
        } else if (ipack == 7) {
            //
            if (isym != 1) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j - kuu; i <= j; i = i + 1) {
                        a[((i - j + kuu + 1) - 1) + (j - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                        if (i < 1) {
                            a[((j - i + 1 + kuu) - 1) + ((i + n) - 1) * lda] = czero;
                        }
                        if (i >= 1 && i != j) {
                            if (isym == 0) {
                                a[((j - i + 1 + kuu) - 1) + (i - 1) * lda] = conj(a[((i - j + kuu + 1) - 1) + (j - 1) * lda]);
                            } else {
                                a[((j - i + 1 + kuu) - 1) + (i - 1) * lda] = a[((i - j + kuu + 1) - 1) + (j - 1) * lda];
                            }
                        }
                    }
                }
            } else if (isym == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j - kuu; i <= j + kll; i = i + 1) {
                        a[((i - j + kuu + 1) - 1) + (j - 1) * lda] = Clatm2(m, n, i, j, kl, ku, idist, iseed, d, igrade, dl, dr, ipvtng, iwork, sparse);
                    }
                }
            }
            //
        }
        //
    }
    //
    //     5)      Scaling the norm
    //
    REAL tempa[1];
    REAL onorm = 0.0;
    if (ipack == 0) {
        onorm = Clange("M", m, n, a, lda, tempa);
    } else if (ipack == 1) {
        onorm = Clansy("M", "U", n, a, lda, tempa);
    } else if (ipack == 2) {
        onorm = Clansy("M", "L", n, a, lda, tempa);
    } else if (ipack == 3) {
        onorm = Clansp("M", "U", n, a, tempa);
    } else if (ipack == 4) {
        onorm = Clansp("M", "L", n, a, tempa);
    } else if (ipack == 5) {
        onorm = Clansb("M", "L", n, kll, a, lda, tempa);
    } else if (ipack == 6) {
        onorm = Clansb("M", "U", n, kuu, a, lda, tempa);
    } else if (ipack == 7) {
        onorm = Clangb("M", n, kll, kuu, a, lda, tempa);
    }
    //
    if (anorm >= zero) {
        //
        if (anorm > zero && onorm == zero) {
            //
            //           Desired scaling impossible
            //
            info = 5;
            return;
            //
        } else if ((anorm > one && onorm < one) || (anorm < one && onorm > one)) {
            //
            //           Scale carefully to avoid over / underflow
            //
            if (ipack <= 2) {
                for (j = 1; j <= n; j = j + 1) {
                    CRscal(m, one / onorm, &a[(j - 1) * lda], 1);
                    CRscal(m, anorm, &a[(j - 1) * lda], 1);
                }
                //
            } else if (ipack == 3 || ipack == 4) {
                //
                CRscal(n * (n + 1) / 2, one / onorm, a, 1);
                CRscal(n * (n + 1) / 2, anorm, a, 1);
                //
            } else if (ipack >= 5) {
                //
                for (j = 1; j <= n; j = j + 1) {
                    CRscal(kll + kuu + 1, one / onorm, &a[(j - 1) * lda], 1);
                    CRscal(kll + kuu + 1, anorm, &a[(j - 1) * lda], 1);
                }
                //
            }
            //
        } else {
            //
            //           Scale straightforwardly
            //
            if (ipack <= 2) {
                for (j = 1; j <= n; j = j + 1) {
                    CRscal(m, anorm / onorm, &a[(j - 1) * lda], 1);
                }
                //
            } else if (ipack == 3 || ipack == 4) {
                //
                CRscal(n * (n + 1) / 2, anorm / onorm, a, 1);
                //
            } else if (ipack >= 5) {
                //
                for (j = 1; j <= n; j = j + 1) {
                    CRscal(kll + kuu + 1, anorm / onorm, &a[(j - 1) * lda], 1);
                }
            }
            //
        }
        //
    }
    //
    //     End of Clatmr
    //
}
